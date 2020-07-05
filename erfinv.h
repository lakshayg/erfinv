#pragma once

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

// Returns a floating point number y such that std::erf(y)
// is close to x. The current implementation is quite accurate
// when x is away from +1.0 and -1.0. As x approaches closer
// to those values, the error in the result increases.
EXTERN_C long double erfinv(long double x);

// Refine the result of erfinv by performing Newton-Raphson
// iteration nr_iter number of times. This method works well
// when the value of x is away from 1.0 and -1.0
EXTERN_C long double erfinv_refine(long double x, int nr_iter);
