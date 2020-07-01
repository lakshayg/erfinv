#include "ppnd.h"
#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"

template <typename T>
class range {
 public:
  range(T value, T eps)
    : value_(value)
    , eps_(eps) {}

  friend bool operator==(T lhs, const range<T>& rhs) {
    return (lhs >= rhs.value_ - rhs.eps_) && (lhs <= rhs.value_ + rhs.eps_);
  }

  friend bool operator==(const range<T>& lhs, T rhs) {
    return (rhs >= lhs.value_ - lhs.eps_) && (rhs <= lhs.value_ + lhs.eps_);
  }

  friend std::ostream& operator<<(std::ostream& os, const range<T>& range) {
    return os << range.value_ << " ± " << range.eps_;
  }

 private:
  T value_, eps_;
};

template <typename T>
range<T> approx(T value, T eps) {
  return range<T>(value, eps);
}

TEST_CASE("PPND7", "Suggested cases") {
  SECTION("float") {
    REQUIRE(ppnd7(0.25f) == approx(-0.6744897501960817f, 1e-7f));
    REQUIRE(ppnd7(1e-3f) == approx(-3.090232306167814f, 1e-6f));
    REQUIRE(ppnd7(1e-20f) == approx(-9.262340089798408f, 1e-8f));
  }

  SECTION("double") {
    REQUIRE(ppnd7(0.25) == approx(-0.6744897501960817, 1e-7));
    REQUIRE(ppnd7(1e-3) == approx(-3.090232306167814, 1e-6));
    REQUIRE(ppnd7(1e-20) == approx(-9.262340089798408, 1e-6));
  }

  SECTION("long double") {
    REQUIRE(ppnd7(0.25l) == approx(-0.6744897501960817l, 1e-7l));
    REQUIRE(ppnd7(1e-3l) == approx(-3.090232306167814l, 1e-6l));
    REQUIRE(ppnd7(1e-20l) == approx(-9.262340089798408l, 1e-6l));
  }
}

TEST_CASE("PPND16", "Suggested cases") {
  SECTION("float") {
    REQUIRE(ppnd16(0.25f) == approx(-0.6744897501960817f, 1e-7f));
    REQUIRE(ppnd16(1e-3f) == approx(-3.090232306167814f, 1e-6f));
    REQUIRE(ppnd16(1e-20f) == approx(-9.262340089798408f, 1e-6f));
  }

  SECTION("double") {
    REQUIRE(ppnd16(0.25) == approx(-0.6744897501960817, 1e-16));
    REQUIRE(ppnd16(1e-3) == approx(-3.090232306167814, 1e-15));
    REQUIRE(ppnd16(1e-20) == approx(-9.262340089798408, 1e-15));
  }

  SECTION("long double") {
    REQUIRE(ppnd16(0.25l) == approx(-0.6744897501960817l, 1e-16l));
    REQUIRE(ppnd16(1e-3l) == approx(-3.090232306167814l, 1e-15l));
    REQUIRE(ppnd16(1e-20l) == approx(-9.262340089798408l, 1e-15l));
  }
}
