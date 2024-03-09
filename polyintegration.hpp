#ifndef HEADERFILE_POLYINTHEADER
#define HEADERFILE_POLYINTHEADER

#include <vector>

namespace poly {
    template <typename T>
    double linearIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata);
    
    template <typename T>
    double linearIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata, T lowLimit, T highLimit);
    
    template <typename T>
    double cubicIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata);
    
    template <typename T>
    double cubicIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata, T lowLimit, T highLimit);
}

#endif