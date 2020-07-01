#pragma once

#include <cassert>
#include <cmath>

namespace detail {
template <typename T>
struct ppnd7_impl {
  static constexpr T a0 = 3.3871327179E0L;
  static constexpr T a1 = 5.0434271938E1L;
  static constexpr T a2 = 1.5929113202E2L;
  static constexpr T a3 = 5.9109374720E1L;

  static constexpr T b1 = 1.7895169469E1L;
  static constexpr T b2 = 7.8757757664E1L;
  static constexpr T b3 = 6.7187563600E1L;

  static constexpr T c0 = 1.4234372777E0L;
  static constexpr T c1 = 2.7568153900E0L;
  static constexpr T c2 = 1.3067284816E0L;
  static constexpr T c3 = 1.7023821103E-1L;

  static constexpr T d1 = 7.3700164250E-1L;
  static constexpr T d2 = 1.2021132975E-1L;

  static constexpr T e0 = 6.6579051150E0L;
  static constexpr T e1 = 3.0812263860E0L;
  static constexpr T e2 = 4.2868294337E-1L;
  static constexpr T e3 = 1.7337203997E-2L;

  static constexpr T f1 = 2.4197894225E-1L;
  static constexpr T f2 = 1.2258202635E-2L;

  static inline T approx1(T r) {
    T x0 = a0 + a1 * r;
    T x1 = a2 + a3 * r;
    T y0 = 1. + b1 * r;
    T y1 = b2 + b3 * r;
    T r2 = r * r;
    return (x0 + x1 * r2) / (y0 + y1 * r2);
  }

  static inline T approx2(T r) {
    T x0 = c0 + c1 * r;
    T x1 = c2 + c3 * r;
    T y0 = d1 + d2 * r;
    return (x0 + x1 * r * r) / (1. + y0 * r);
  }

  static inline T approx3(T r) {
    T x0 = e0 + e1 * r;
    T x1 = e2 + e3 * r;
    T y0 = f1 + f2 * r;
    return (x0 + x1 * r * r) / (1. + y0 * r);
  }
};

template <typename T>
struct ppnd16_impl {
  static constexpr T a0 = 3.3871328727963666080E0L;
  static constexpr T a1 = 1.3314166789178437745E2L;
  static constexpr T a2 = 1.9715909503065514427E3L;
  static constexpr T a3 = 1.3731693765509461125E4L;
  static constexpr T a4 = 4.5921953931549871457E4L;
  static constexpr T a5 = 6.7265770927008700853E4L;
  static constexpr T a6 = 3.3430575583588128105E4L;
  static constexpr T a7 = 2.5090809287301226727E3L;

  static constexpr T b1 = 4.2313330701600911252E1L;
  static constexpr T b2 = 6.8718700749205790830E2L;
  static constexpr T b3 = 5.3941960214247511077E3L;
  static constexpr T b4 = 2.1213794301586595867E4L;
  static constexpr T b5 = 3.9307895800092710610E4L;
  static constexpr T b6 = 2.8729085735721942674E4L;
  static constexpr T b7 = 5.2264952788528545610E3L;

  static constexpr T c0 = 1.42343711074968357734E0L;
  static constexpr T c1 = 4.63033784615654529590E0L;
  static constexpr T c2 = 5.76949722146069140550E0L;
  static constexpr T c3 = 3.64784832476320460504E0L;
  static constexpr T c4 = 1.27045825245236838258E0L;
  static constexpr T c5 = 2.41780725177450611770E-1L;
  static constexpr T c6 = 2.27238449892691845833E-2L;
  static constexpr T c7 = 7.74545014278341407640E-4L;

  static constexpr T d1 = 2.05319162663775882187E0L;
  static constexpr T d2 = 1.67638483018380384940E0L;
  static constexpr T d3 = 6.89767334985100004550E-1L;
  static constexpr T d4 = 1.48103976427480074590E-1L;
  static constexpr T d5 = 1.51986665636164571966E-2L;
  static constexpr T d6 = 5.47593808499534494600E-4L;
  static constexpr T d7 = 1.05075007164441684324E-9L;

  static constexpr T e0 = 6.65790464350110377720E0L;
  static constexpr T e1 = 5.46378491116411436990E0L;
  static constexpr T e2 = 1.78482653991729133580E0L;
  static constexpr T e3 = 2.96560571828504891230E-1L;
  static constexpr T e4 = 2.65321895265761230930E-2L;
  static constexpr T e5 = 1.24266094738807843860E-3L;
  static constexpr T e6 = 2.71155556874348757815E-5L;
  static constexpr T e7 = 2.01033439929228813265E-7L;

  static constexpr T f1 = 5.99832206555887937690E-1L;
  static constexpr T f2 = 1.36929880922735805310E-1L;
  static constexpr T f3 = 1.48753612908506148525E-2L;
  static constexpr T f4 = 7.86869131145613259100E-4L;
  static constexpr T f5 = 1.84631831751005468180E-5L;
  static constexpr T f6 = 1.42151175831644588870E-7L;
  static constexpr T f7 = 2.04426310338993978564E-15L;

  static inline T approx1(T r) {
    T x00 = a0 + a1 * r;
    T x01 = a2 + a3 * r;
    T x02 = a4 + a5 * r;
    T x03 = a6 + a7 * r;
    T y00 = 1. + b1 * r;
    T y01 = b2 + b3 * r;
    T y02 = b4 + b5 * r;
    T y03 = b6 + b7 * r;

    T r2 = r * r;
    T x10 = x00 + x01 * r2;
    T x11 = x02 + x03 * r2;
    T y10 = y00 + y01 * r2;
    T y11 = y02 + y03 * r2;

    T r4 = r2 * r2;
    T x20 = x10 + x11 * r4;
    T y20 = y10 + y11 * r4;

    return x20 / y20;
  }

  static inline T approx2(T r) {
    T x00 = c0 + c1 * r;
    T x01 = c2 + c3 * r;
    T x02 = c4 + c5 * r;
    T x03 = c6 + c7 * r;
    T y00 = 1. + d1 * r;
    T y01 = d2 + d3 * r;
    T y02 = d4 + d5 * r;
    T y03 = d6 + d7 * r;

    T r2 = r * r;
    T x10 = x00 + x01 * r2;
    T x11 = x02 + x03 * r2;
    T y10 = y00 + y01 * r2;
    T y11 = y02 + y03 * r2;

    T r4 = r2 * r2;
    T x20 = x10 + x11 * r4;
    T y20 = y10 + y11 * r4;

    return x20 / y20;
  }

  static inline T approx3(T r) {
    T x00 = e0 + e1 * r;
    T x01 = e2 + e3 * r;
    T x02 = e4 + e5 * r;
    T x03 = e6 + e7 * r;
    T y00 = 1. + f1 * r;
    T y01 = f2 + f3 * r;
    T y02 = f4 + f5 * r;
    T y03 = f6 + f7 * r;

    T r2 = r * r;
    T x10 = x00 + x01 * r2;
    T x11 = x02 + x03 * r2;
    T y10 = y00 + y01 * r2;
    T y11 = y02 + y03 * r2;

    T r4 = r2 * r2;
    T x20 = x10 + x11 * r4;
    T y20 = y10 + y11 * r4;

    return x20 / y20;
  }
};

template <typename T, typename impl>
inline T ppnd_main(T p) {
  assert(p >= 0.0 && p <= 1.0);
  T q = p - 0.5;

  if (std::fabs(q) <= 0.425) {
    T r = 0.180625 - q * q;
    return q * impl::approx1(r);
  }

  T u = q < 0.0 ? p : (1.0 - p);
  T r = std::sqrt(-std::log(u));
  T approx = r > 5.0 ? impl::approx3(r - 5.0) : impl::approx2(r - 1.6);
  return std::copysign(approx, q);
}
} // namespace detail

template <typename T>
inline T ppnd7(T p) {
  return detail::ppnd_main<T, detail::ppnd7_impl<T>>(p);
}

template <typename T>
inline T ppnd16(T p) {
  return detail::ppnd_main<T, detail::ppnd16_impl<T>>(p);
}
