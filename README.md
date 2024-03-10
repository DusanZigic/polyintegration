# POLYINTEGRATE  
Set of functions that integrate interpolated polynomial for given data. Interpolated polynomial can be linear or cubic.  Functions are templated to account for float, double and long double data.  

Integration can be done on the entire range of the data:  
<pre>
 T linearIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata)
 T cubicIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata)
</pre>
or, integration limits can be specified:  
<pre>
 T linearIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata, T lowLimit, T highLimit)
 T cubicIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata, T lowLimit, T highLimit)
</pre>

All 4 functions are embeded in *poly* namespace.  

Sample usage to integrate *f(x)=2x* from 1.5 to 3.5 would be:
```
std::vector<double> x{1.0, 2.0, 3.0, 4.0}, f{2.0, 4.0, 6.0, 8.0};
double result = poly::linearIntegrate(x, f, 1.5, 3.5);
```
