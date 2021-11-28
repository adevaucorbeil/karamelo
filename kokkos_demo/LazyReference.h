#pragma once

template<typename T>
class LazyReference {
    T *source;
    T local;

public:
    LazyReference() = default;

    LazyReference(T *source):
        source(source),
        local(*source)
    {}
    
    ~LazyReference() {
        if (*source != local) {
            *source = local;
        }
    }

    operator const T &() const {
      return local;
    }

    operator T &() {
        return local;
    }

    T &
    operator=(const T &value)
    {
      return local = value;
    }

    T &
    operator/=(const T &value)
    {
      return local /= value;
    }

    T &
    operator*=(const T &value)
    {
      return local *= value;
    }
};