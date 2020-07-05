Implementation of the inverse error function based on the rational
approximation of percentage points of normal distribution available from
https://www.jstor.org/stable/2347330.

                1           /  x + 1  \
erfinv(x) = --------- ppnd |  -------  |
             sqrt(2)        \    2    /

The code has been tested on an x86_64 machine with Intel Core i7, the
tests provided in this repository might not pass for different hardware
configuration due to difference in floating point operations.

golang's math library uses the same implementation for erfinv
