// A small stand-in for geogram/basic/attributes.h.
//
// Reimplemented rather than copied: an attribute store is addressed by element index, so it has no
// tie-breaking of its own. Only the slice the mesh container and fTetWild use is here -- named
// attributes of a scalar type, bound or created on first use, resized and permuted with their
// elements.

#pragma once

#include "geo_basic.h"

#include <map>
#include <memory>
#include <string>

namespace floatTetWild {
namespace geo {

    /**
     * \brief Base class for a named array of per-element values.
     */
    class AttributeStoreBase {
    public:
        virtual ~AttributeStoreBase() {
        }

        virtual void resize(index_t nb) = 0;

        /**
         * \brief Reorders the values so that value \p k becomes the value that was at
         *  \p permutation[k].
         */
        virtual void apply_permutation(const vector<index_t>& permutation) = 0;
    };

    /**
     * \brief A named array of per-element values of type \p T.
     */
    template <class T>
    class AttributeStore : public AttributeStoreBase {
    public:
        void resize(index_t nb) override {
            data_.resize(nb, T());
        }

        void apply_permutation(const vector<index_t>& permutation) override {
            geo_debug_assert(permutation.size() <= data_.size());
            std::vector<T> reordered(data_.size());
            for(index_t k = 0; k < permutation.size(); ++k) {
                reordered[k] = data_[permutation[k]];
            }
            data_.swap(reordered);
        }

        T& operator[] (index_t i) {
            geo_debug_assert(i < data_.size());
            return data_[i];
        }

        const T& operator[] (index_t i) const {
            geo_debug_assert(i < data_.size());
            return data_[i];
        }

    private:
        std::vector<T> data_;
    };

    /**
     * \brief The attributes attached to one kind of mesh element.
     */
    class AttributesManager {
    public:
        index_t size() const {
            return size_;
        }

        void resize(index_t nb) {
            for(auto& entry : stores_) {
                entry.second->resize(nb);
            }
            size_ = nb;
        }

        void apply_permutation(const vector<index_t>& permutation) {
            for(auto& entry : stores_) {
                entry.second->apply_permutation(permutation);
            }
        }

        void clear(bool keep_attributes) {
            if(keep_attributes) {
                resize(0);
            } else {
                stores_.clear();
                size_ = 0;
            }
        }

        /**
         * \brief Finds the store named \p name, creating it if it does not exist yet.
         */
        template <class T>
        AttributeStore<T>* find_or_create(const std::string& name) {
            auto it = stores_.find(name);
            if(it != stores_.end()) {
                return dynamic_cast<AttributeStore<T>*>(it->second.get());
            }
            AttributeStore<T>* store = new AttributeStore<T>;
            store->resize(size_);
            stores_[name].reset(store);
            return store;
        }

    private:
        std::map<std::string, std::unique_ptr<AttributeStoreBase>> stores_;
        index_t size_ = 0;
    };

    /**
     * \brief A reference to a named attribute, bound to it or creating it on construction.
     */
    template <class T>
    class Attribute {
    public:
        Attribute(AttributesManager& manager, const std::string& name) :
            store_(manager.find_or_create<T>(name)) {
            geo_assert(store_ != nullptr);
        }

        T& operator[] (index_t i) {
            return (*store_)[i];
        }

        const T& operator[] (index_t i) const {
            return (*store_)[i];
        }

        T& operator[] (int i) {
            return (*store_)[index_t(i)];
        }

        const T& operator[] (int i) const {
            return (*store_)[index_t(i)];
        }

    private:
        AttributeStore<T>* store_;
    };
}
}
