#pragma once

template<typename T>
class LazyReference {
    T *source;
    T local;

public:
    LazyReference(T *source):
        source(source),
        local(*source)
    {}
    
    ~LazyReference() {
        if (*source != local) {
            *source = local;
        }
    }

    operator T &() const {
        return local;
    }
};