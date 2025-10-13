//
// SPDX-FileCopyrightText: 2017-2025 Lakshay Garg
// SPDX-License-Identifier: MIT
//

#include "erfinv.h"
#include <math.h>

EXTERN_C long double erfinv(long double x) {

  if (x < -1 || x > 1) {
    return NAN;
  } else if (x == 1.0) {
    return INFINITY;
  } else if (x == -1.0) {
    return -INFINITY;
  }

  const long double LN2 = 6.931471805599453094172321214581e-1L;

  const long double A0 = 1.1975323115670912564578e0L;
  const long double A1 = 4.7072688112383978012285e1L;
  const long double A2 = 6.9706266534389598238465e2L;
  const long double A3 = 4.8548868893843886794648e3L;
  const long double A4 = 1.6235862515167575384252e4L;
  const long double A5 = 2.3782041382114385731252e4L;
  const long double A6 = 1.1819493347062294404278e4L;
  const long double A7 = 8.8709406962545514830200e2L;

  const long double B0 = 1.0000000000000000000e0L;
  const long double B1 = 4.2313330701600911252e1L;
  const long double B2 = 6.8718700749205790830e2L;
  const long double B3 = 5.3941960214247511077e3L;
  const long double B4 = 2.1213794301586595867e4L;
  const long double B5 = 3.9307895800092710610e4L;
  const long double B6 = 2.8729085735721942674e4L;
  const long double B7 = 5.2264952788528545610e3L;

  const long double C0 = 1.42343711074968357734e0L;
  const long double C1 = 4.63033784615654529590e0L;
  const long double C2 = 5.76949722146069140550e0L;
  const long double C3 = 3.64784832476320460504e0L;
  const long double C4 = 1.27045825245236838258e0L;
  const long double C5 = 2.41780725177450611770e-1L;
  const long double C6 = 2.27238449892691845833e-2L;
  const long double C7 = 7.74545014278341407640e-4L;

  const long double D0 = 1.4142135623730950488016887e0L;
  const long double D1 = 2.9036514445419946173133295e0L;
  const long double D2 = 2.3707661626024532365971225e0L;
  const long double D3 = 9.7547832001787427186894837e-1L;
  const long double D4 = 2.0945065210512749128288442e-1L;
  const long double D5 = 2.1494160384252876777097297e-2L;
  const long double D6 = 7.7441459065157709165577218e-4L;
  const long double D7 = 1.4859850019840355905497876e-9L;

  const long double E0 = 6.65790464350110377720e0L;
  const long double E1 = 5.46378491116411436990e0L;
  const long double E2 = 1.78482653991729133580e0L;
  const long double E3 = 2.96560571828504891230e-1L;
  const long double E4 = 2.65321895265761230930e-2L;
  const long double E5 = 1.24266094738807843860e-3L;
  const long double E6 = 2.71155556874348757815e-5L;
  const long double E7 = 2.01033439929228813265e-7L;

  const long double F0 = 1.414213562373095048801689e0L;
  const long double F1 = 8.482908416595164588112026e-1L;
  const long double F2 = 1.936480946950659106176712e-1L;
  const long double F3 = 2.103693768272068968719679e-2L;
  const long double F4 = 1.112800997078859844711555e-3L;
  const long double F5 = 2.611088405080593625138020e-5L;
  const long double F6 = 2.010321207683943062279931e-7L;
  const long double F7 = 2.891024605872965461538222e-15L;

  long double abs_x = fabsl(x);
  long double r = 0.0L;
  long double num = 0.0L;
  long double den = 0.0L;

  if (abs_x <= 0.85L) {
    r = 0.180625L - 0.25L * x * x;

    num = fmal(num, r, A7);
    num = fmal(num, r, A6);
    num = fmal(num, r, A5);
    num = fmal(num, r, A4);
    num = fmal(num, r, A3);
    num = fmal(num, r, A2);
    num = fmal(num, r, A1);
    num = fmal(num, r, A0);

    den = fmal(den, r, B7);
    den = fmal(den, r, B6);
    den = fmal(den, r, B5);
    den = fmal(den, r, B4);
    den = fmal(den, r, B3);
    den = fmal(den, r, B2);
    den = fmal(den, r, B1);
    den = fmal(den, r, B0);

    return x * num / den;
  }

  r = sqrtl(LN2 - logl(1.0L - abs_x));

  if (r <= 5.0L) {
    r = r - 1.6L;

    num = fmal(num, r, C7);
    num = fmal(num, r, C6);
    num = fmal(num, r, C5);
    num = fmal(num, r, C4);
    num = fmal(num, r, C3);
    num = fmal(num, r, C2);
    num = fmal(num, r, C1);
    num = fmal(num, r, C0);

    den = fmal(den, r, D7);
    den = fmal(den, r, D6);
    den = fmal(den, r, D5);
    den = fmal(den, r, D4);
    den = fmal(den, r, D3);
    den = fmal(den, r, D2);
    den = fmal(den, r, D1);
    den = fmal(den, r, D0);
  } else {
    r = r - 5.0L;

    num = fmal(num, r, E7);
    num = fmal(num, r, E6);
    num = fmal(num, r, E5);
    num = fmal(num, r, E4);
    num = fmal(num, r, E3);
    num = fmal(num, r, E2);
    num = fmal(num, r, E1);
    num = fmal(num, r, E0);

    den = fmal(den, r, F7);
    den = fmal(den, r, F6);
    den = fmal(den, r, F5);
    den = fmal(den, r, F4);
    den = fmal(den, r, F3);
    den = fmal(den, r, F2);
    den = fmal(den, r, F1);
    den = fmal(den, r, F0);
  }

  return copysignl(num / den, x);
}

EXTERN_C long double erfinv_refine(long double x, int nr_iter) {
  const long double k = 0.8862269254527580136490837416706L; // 0.5 * sqrt(pi)
  long double y = erfinv(x);
  while (nr_iter-- > 0) {
    y -= k * (erfl(y) - x) / expl(-y * y);
  }
  return y;
}
