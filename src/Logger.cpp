// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Jeremie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Jeremie Dumas on 09/04/18.
//

#include "Logger.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace floatTetWild {
namespace {

// spdlog's names and console colours, so a run looks the same as it did.
const char* level_name(Logger::Level level)
{
    switch (level) {
    case Logger::Trace:
        return "trace";
    case Logger::Debug:
        return "debug";
    case Logger::Info:
        return "info";
    case Logger::Warn:
        return "warning";
    case Logger::Error:
        return "error";
    case Logger::Critical:
        return "critical";
    case Logger::Off:
    default:
        return "off";
    }
}

const char* level_colour(Logger::Level level)
{
    switch (level) {
    case Logger::Trace:
        return "\033[37m";
    case Logger::Debug:
        return "\033[36m";
    case Logger::Info:
        return "\033[32m";
    case Logger::Warn:
        return "\033[33m\033[1m";
    case Logger::Error:
        return "\033[31m\033[1m";
    case Logger::Critical:
        return "\033[1m\033[41m";
    case Logger::Off:
    default:
        return "\033[m";
    }
}

// Colour only when stdout is a terminal, which is what spdlog's console sink did. It keeps
// escape codes out of redirected output.
bool stdout_is_terminal()
{
#ifdef _WIN32
    return _isatty(_fileno(stdout)) != 0;
#else
    return isatty(fileno(stdout)) != 0;
#endif
}

std::string timestamp()
{
    using namespace std::chrono;
    const auto now  = system_clock::now();
    const auto secs = system_clock::to_time_t(now);
    const auto ms   = duration_cast<milliseconds>(now.time_since_epoch()) % seconds(1);

    std::tm broken_down {};
#ifdef _WIN32
    localtime_s(&broken_down, &secs);
#else
    localtime_r(&secs, &broken_down);
#endif
    char date[32];
    std::strftime(date, sizeof(date), "%Y-%m-%d %H:%M:%S", &broken_down);

    std::ostringstream out;
    out << date << '.' << std::setfill('0') << std::setw(3) << ms.count();
    return out.str();
}

}  // namespace

Logger& logger()
{
    static Logger instance;
    return instance;
}

void Logger::init(bool use_cout, const std::string& filename, bool truncate)
{
    Logger&                     instance = logger();
    std::lock_guard<std::mutex> lock(instance.mutex_);

    instance.use_cout_ = use_cout;
    if (instance.file_.is_open())
        instance.file_.close();
    if (!filename.empty())
        instance.file_.open(filename, truncate ? std::ios::trunc : std::ios::app);
}

void Logger::set_level(int level)
{
    std::lock_guard<std::mutex> lock(mutex_);
    level_ = level < Trace ? int(Trace) : (level > Off ? int(Off) : level);
}

void Logger::write(Level level, const std::string& message)
{
    // Only the level name is coloured, matching spdlog's default pattern.
    static const bool colour = stdout_is_terminal();

    const std::string stamp  = timestamp();
    const std::string prefix = "[" + stamp + "] [float-tetwild] [";

    // Flushed per line. spdlog logged from a background thread and flushed on a timer, so its
    // lines could land out of order against the plain cout printing the rest of the code does.
    // Writing them on the calling thread keeps the two in the order they were issued.
    std::lock_guard<std::mutex> lock(mutex_);
    if (use_cout_) {
        if (colour)
            std::cout << prefix << level_colour(level) << level_name(level) << "\033[m] " << message
                      << std::endl;
        else
            std::cout << prefix << level_name(level) << "] " << message << std::endl;
    }
    if (file_.is_open())
        file_ << prefix << level_name(level) << "] " << message << std::endl;
}

}  // namespace floatTetWild
