#include "erfinv.h"
#include <assert.h>



////////////////////////////////////////////////////////////////////////////////
// data for testing (computed using wolfram alpha)

// A vector containing (x, expected erfinv, allowed error)
//using TABLE = vector<tuple<long double, long double, long double>>;

// Nearly uniform samples of x from [-1, 1]
const long double uniform_data[23][3] = {
  {-0.9999901234567890L, -3.12531209890571959826L, 1e-15L}, // 0
  {-0.9090817727272728L, -1.19541621600439254353L, 1e-16L}, // 1
  {-0.8181735454545455L, -0.94409535430005211436L, 1e-16L}, // 2
  {-0.7272653181818183L, -0.77554525438923442604L, 1e-16L}, // 3
  {-0.6363570909090909L, -0.64236795545430416665L, 1e-16L}, // 4
  {-0.5454488636363637L, -0.52880922364399673143L, 1e-16L}, // 5
  {-0.4545406363636364L, -0.42750127212741209472L, 1e-16L}, // 6
  {-0.3636324090909091L, -0.33430847448420056223L, 1e-16L}, // 7
  {-0.2727241818181818L, -0.24660460627491576384L, 1e-16L}, // 8
  {-0.1818159545454545L, -0.16255059170754568956L, 1e-17L}, // 9
  {-0.0909077272727272L, -0.08073997953126072744L, 1e-17L}, // 10
  { 0.00000050123456789L, 0.00000044420757003182L, 1e-19L}, // 11
  { 0.09090872727272725L, 0.08074087155438624655L, 1e-17L}, // 12
  { 0.18181695454545455L, 0.16255150166321581410L, 1e-17L}, // 13
  { 0.27272518181818184L, 0.24660554806942288638L, 1e-17L}, // 14
  { 0.36363340909090913L, 0.33430946550500629626L, 1e-16L}, // 15
  { 0.45454163636363636L, 0.42750233606374840004L, 1e-16L}, // 16
  { 0.54544986363636364L, 0.52881039581498761461L, 1e-16L}, // 17
  { 0.63635809090909090L, 0.64236929436345064007L, 1e-16L}, // 18
  { 0.72726631818181831L, 0.77554687157823653280L, 1e-16L}, // 19
  { 0.81817454545454562L, 0.94409751522390649338L, 1e-16L}, // 20
  { 0.90908277272727291L, 1.19541991566879908234L, 1e-16L}, // 21
  { 0.99999100123456789L, 3.13950543736986052401L, 1e-15L}, // 22
};

// Values of x very close to the asymptotes at x=+1 and x=-1
const long double  close_to_1[6][3] = {
  {-0.99999999999L, -4.8129240673658227430L, 1e-09L},
  {-0.99999999990L, -4.5728249673894852787L, 1e-10L},
  {-0.99999999900L, -4.3200053849134452862L, 1e-11L},
  { 0.99999999900L,  4.3200053849134452862L, 1e-11L},
  { 0.99999999990L,  4.5728249673894852787L, 1e-10L},
  { 0.99999999999L,  4.8129240673658227430L, 1e-09L},
}; // size = 6

// The implementation uses functions to approximate erfinv in
// various intervals. This table contains the points where two
// different functions come together
const long double  joining_points[6][3] = {
  { 0.85000000000L,  1.0179024648320276436L, 1e-16},
  {-0.85000000000L, -1.0179024648320276436L, 1e-16},
  // |x| is smaller than 1 - 2 * exp(-25)
  { 0.99999999997L,  4.6998353461721558307L, 1e-10},
  {-0.99999999997L, -4.6998353461721558307L, 1e-10},
  // |x| is slightly larger than 1 - 2 * exp(-25)
  { 0.99999999998L,  4.7418744480446202994L, 1e-09},
  {-0.99999999998L, -4.7418744480446202994L, 1e-09},
}; // size = 6

////////////////////////////////////////////////////////////////////////////////
// Unit tests

#define ASSERT_INF(x, sign)                       \
  do {                                            \
    assert(sign != 0.0);                            \
    long double y = (x);                             \
    if (sign > 0) {                               \
      assert(isinf(y) && !signbit(y));            \
    } else if (sign < 0) {                        \
      assert(isinf(y) && signbit(y));             \
    }                                             \
  } while (0)

#define ASSERT_NAN(x) assert(isnan(x))

void test_out_of_domain_values() {
  ASSERT_NAN(erfinv( (long double) -2.0));
  ASSERT_INF(erfinv( (long double) -1.0), (long double) -1.0);
  ASSERT_INF(erfinv( (long double)  1.0), (long double)  1.0);
  ASSERT_NAN(erfinv( (long double)  2.0));
}

#define ASSERT_ERFINV(table_entry) {              \
  long double x, y, eps;                          \
  x=table_entry[0];                               \
  y=table_entry[1];                               \
  eps=table_entry[2];                             \
  assert(fabsl(erfinv(x) - y) < eps);             \
}

void test_values_close_to_asymptotes() {
  ASSERT_ERFINV(close_to_1[0]);
  ASSERT_ERFINV(close_to_1[1]);
  ASSERT_ERFINV(close_to_1[2]);
  ASSERT_ERFINV(close_to_1[3]);
  ASSERT_ERFINV(close_to_1[4]);
  ASSERT_ERFINV(close_to_1[5]);
}

void test_piecewise_boundaries() {
  ASSERT_ERFINV(joining_points[0]);
  ASSERT_ERFINV(joining_points[1]);
  ASSERT_ERFINV(joining_points[2]);
  ASSERT_ERFINV(joining_points[3]);
  ASSERT_ERFINV(joining_points[4]);
  ASSERT_ERFINV(joining_points[5]);
}

void test_uniform_samples() {
  ASSERT_ERFINV(uniform_data[0]);
  ASSERT_ERFINV(uniform_data[1]);
  ASSERT_ERFINV(uniform_data[2]);
  ASSERT_ERFINV(uniform_data[3]);
  ASSERT_ERFINV(uniform_data[4]);
  ASSERT_ERFINV(uniform_data[5]);
  ASSERT_ERFINV(uniform_data[6]);
  ASSERT_ERFINV(uniform_data[7]);
  ASSERT_ERFINV(uniform_data[8]);
  ASSERT_ERFINV(uniform_data[9]);
  ASSERT_ERFINV(uniform_data[10]);
  ASSERT_ERFINV(uniform_data[11]);
  ASSERT_ERFINV(uniform_data[12]);
  ASSERT_ERFINV(uniform_data[13]);
  ASSERT_ERFINV(uniform_data[14]);
  ASSERT_ERFINV(uniform_data[15]);
  ASSERT_ERFINV(uniform_data[16]);
  ASSERT_ERFINV(uniform_data[17]);
  ASSERT_ERFINV(uniform_data[18]);
  ASSERT_ERFINV(uniform_data[19]);
  ASSERT_ERFINV(uniform_data[20]);
  ASSERT_ERFINV(uniform_data[21]);
  ASSERT_ERFINV(uniform_data[22]);
}

#define DISABLE_ASSERT_ERFINV_WITH_NR(x) (void)0
#define ASSERT_ERFINV_WITH_NR(table_entry) {      \
  long double x, y, eps;                          \
  x  =table_entry[0];                             \
  y  =table_entry[1];                             \
  eps=table_entry[2];                             \
  assert(fabsl(erfinv_iter(x, 100000) - y) < 0.1 * eps);    \
}

void test_erfinv_with_newton_raphson() {
  // removing points at the extremeties as NR is unable
  // to compute a better estimate for those points
  DISABLE_ASSERT_ERFINV_WITH_NR(uniform_data[0]);
  ASSERT_ERFINV_WITH_NR(uniform_data[1]);
  ASSERT_ERFINV_WITH_NR(uniform_data[2]);
  ASSERT_ERFINV_WITH_NR(uniform_data[3]);
  ASSERT_ERFINV_WITH_NR(uniform_data[4]);
  ASSERT_ERFINV_WITH_NR(uniform_data[5]);
  ASSERT_ERFINV_WITH_NR(uniform_data[6]);
  ASSERT_ERFINV_WITH_NR(uniform_data[7]);
  ASSERT_ERFINV_WITH_NR(uniform_data[8]);
  ASSERT_ERFINV_WITH_NR(uniform_data[9]);
  ASSERT_ERFINV_WITH_NR(uniform_data[10]);
  ASSERT_ERFINV_WITH_NR(uniform_data[11]);
  ASSERT_ERFINV_WITH_NR(uniform_data[12]);
  ASSERT_ERFINV_WITH_NR(uniform_data[13]);
  ASSERT_ERFINV_WITH_NR(uniform_data[14]);
  ASSERT_ERFINV_WITH_NR(uniform_data[15]);
  ASSERT_ERFINV_WITH_NR(uniform_data[16]);
  ASSERT_ERFINV_WITH_NR(uniform_data[17]);
  ASSERT_ERFINV_WITH_NR(uniform_data[18]);
  ASSERT_ERFINV_WITH_NR(uniform_data[19]);
  ASSERT_ERFINV_WITH_NR(uniform_data[20]);
  ASSERT_ERFINV_WITH_NR(uniform_data[21]);
  DISABLE_ASSERT_ERFINV_WITH_NR(uniform_data[22]);
}

int main() {
printf("Testing out_of_domain_values\n");
  test_out_of_domain_values();
printf("Testing close_to_asymptotes\n");
  test_values_close_to_asymptotes();
printf("Testing piecewise_boundaries\n");
  test_piecewise_boundaries();
printf("Testing uniform_samples\n");
  test_uniform_samples();
printf("Testing erfinv_with_newton_raphson\n");
  test_erfinv_with_newton_raphson();
  return 0;
}
