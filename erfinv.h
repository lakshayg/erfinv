#ifndef LAKSHAYG_ERFINV_H
#define LAKSHAYG_ERFINV_H

// Returns a floating point number y such that std::erf(y)
// is close to x. The current implementation is quite accurate
// when x is away from +1.0 and -1.0. As x approaches closer
// to those values, the error in the result increases.
long double erfinv(long double x);

// Refine the result of erfinv by performing Newton-Raphson
// iteration nr_iter number of times. This method works well
// when the value of x is away from 1.0 and -1.0
long double erfinv(long double x, int nr_iter);

#endif
