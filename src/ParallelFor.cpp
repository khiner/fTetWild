// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/ParallelFor.hpp>

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <thread>

#ifndef _WIN32
#include <pthread.h>
#endif

namespace floatTetWild {
namespace {

// TBB ran its workers on 64 MB stacks and this code still wants them: the exact predicates
// allocate expansions with alloca, and a default pthread stack is 512 KB on macOS.
constexpr size_t WorkerStackSize = 64 * 1024 * 1024;

// Several chunks per thread, so one slow chunk cannot leave the other threads idle. Which
// indices land in which chunk depends only on the range and the thread count, never on
// scheduling.
constexpr size_t ChunksPerThread = 4;

// Set while a thread is running loop bodies. A parallel_for reached from inside one runs
// inline instead of re-entering the pool, which would deadlock on the completion count.
thread_local bool in_parallel_region = false;

struct Range
{
    size_t lo, hi;
};

using Body = std::function<void(size_t, size_t, size_t)>;

class ThreadPool
{
  public:
    explicit ThreadPool(unsigned int parallelism)
        : n_workers_(parallelism > 1 ? parallelism - 1 : 0)
    {
        for (unsigned int i = 0; i < n_workers_; ++i)
            spawn_worker();
    }

    ~ThreadPool()
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        work_ready_.notify_all();
        join_workers();
    }

    // The calling thread runs chunks too, so it counts towards the parallelism.
    size_t parallelism() const { return n_workers_ + 1; }

    void run(const Body& body, const std::vector<Range>& chunks)
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            body_    = &body;
            chunks_  = &chunks;
            next_    = 0;
            running_ = n_workers_;
            ++generation_;
        }
        work_ready_.notify_all();

        claim_chunks();

        std::unique_lock<std::mutex> lock(mutex_);
        work_done_.wait(lock, [this] { return running_ == 0; });
        body_   = nullptr;
        chunks_ = nullptr;
    }

    void worker_loop()
    {
        in_parallel_region = true;
        uint64_t seen      = 0;
        for (;;) {
            {
                std::unique_lock<std::mutex> lock(mutex_);
                work_ready_.wait(lock, [this, seen] { return stopping_ || generation_ != seen; });
                if (stopping_)
                    return;
                seen = generation_;
            }

            claim_chunks();

            std::lock_guard<std::mutex> lock(mutex_);
            if (--running_ == 0)
                work_done_.notify_one();
        }
    }

  private:
    // Runs chunks until none are left. Every thread in a run does this, including the caller.
    void claim_chunks()
    {
        const bool was_in_region = in_parallel_region;
        in_parallel_region       = true;
        for (;;) {
            size_t index;
            Range  range;
            {
                std::lock_guard<std::mutex> lock(mutex_);
                if (chunks_ == nullptr || next_ >= chunks_->size())
                    break;
                index = next_++;
                range = (*chunks_)[index];
            }
            (*body_)(range.lo, range.hi, index);
        }
        in_parallel_region = was_in_region;
    }

    void spawn_worker();
    void join_workers();

    const unsigned int        n_workers_;
    std::mutex                mutex_;
    std::condition_variable   work_ready_, work_done_;
    const Body*               body_       = nullptr;
    const std::vector<Range>* chunks_     = nullptr;
    size_t                    next_       = 0;
    unsigned int              running_    = 0;
    uint64_t                  generation_ = 0;
    bool                      stopping_   = false;

#ifndef _WIN32
    std::vector<pthread_t> workers_;
#else
    std::vector<std::thread> workers_;
#endif
};

#ifndef _WIN32
void* worker_entry(void* self)
{
    static_cast<ThreadPool*>(self)->worker_loop();
    return nullptr;
}

void ThreadPool::spawn_worker()
{
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, WorkerStackSize);
    pthread_t handle;
    if (pthread_create(&handle, &attr, worker_entry, this) == 0)
        workers_.push_back(handle);
    pthread_attr_destroy(&attr);
}

void ThreadPool::join_workers()
{
    for (pthread_t handle : workers_)
        pthread_join(handle, nullptr);
    workers_.clear();
}
#else
// std::thread gives no control over stack size. On Windows the reserve in the executable
// header applies to every thread instead, which CMakeLists sets to match WorkerStackSize.
void ThreadPool::spawn_worker()
{
    workers_.emplace_back([this] { worker_loop(); });
}

void ThreadPool::join_workers()
{
    for (std::thread& worker : workers_)
        if (worker.joinable())
            worker.join();
    workers_.clear();
}
#endif

std::mutex                  pool_mutex;
std::unique_ptr<ThreadPool> pool;
unsigned int                requested_threads = 0;  // 0 asks for one thread per core

ThreadPool& get_pool()
{
    std::lock_guard<std::mutex> lock(pool_mutex);
    if (!pool) {
        const unsigned int n = requested_threads != 0
                                 ? requested_threads
                                 : std::max(1u, std::thread::hardware_concurrency());
        pool.reset(new ThreadPool(n));
    }
    return *pool;
}

std::vector<Range> make_chunks(size_t begin, size_t end, size_t k)
{
    const size_t       n = end - begin;
    std::vector<Range> chunks;
    chunks.reserve(k);
    for (size_t c = 0; c < k; ++c)
        chunks.push_back({begin + (n * c) / k, begin + (n * (c + 1)) / k});
    return chunks;
}

}  // namespace

void set_num_threads(unsigned int n)
{
    std::lock_guard<std::mutex> lock(pool_mutex);
    requested_threads = std::max(1u, n);
    pool.reset();
}

unsigned int get_num_threads()
{ return static_cast<unsigned int>(get_pool().parallelism()); }

size_t detail::chunk_count(size_t begin, size_t end)
{
    if (end <= begin)
        return 0;
    const size_t n = end - begin;
    if (in_parallel_region)
        return 1;
    const size_t threads = get_pool().parallelism();
    if (threads <= 1)
        return 1;
    return std::min(n, threads * ChunksPerThread);
}

void detail::parallel_ranges(size_t begin, size_t end, const Body& body)
{
    const size_t k = chunk_count(begin, end);
    if (k == 0)
        return;
    const std::vector<Range> chunks = make_chunks(begin, end, k);
    if (k == 1) {
        const bool was_in_region = in_parallel_region;
        in_parallel_region       = true;
        body(chunks[0].lo, chunks[0].hi, 0);
        in_parallel_region = was_in_region;
        return;
    }
    get_pool().run(body, chunks);
}

}  // namespace floatTetWild
