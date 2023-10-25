#ifndef REAL_HPP
#define REAL_HPP

#include <cmath>
#include <compare>
#include <concepts>

namespace std {
template <typename T>
concept arithmetic = is_arithmetic<T>::value;
}

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#else
#endif

class unreal_t {
public:
  using num_t = double;
private:
  num_t m_x = 0.;

  static const constexpr num_t Epsilon = 1. / (1l<<40);

public:

  constexpr unreal_t(num_t x = 0.) : m_x(x) {}

  constexpr unreal_t(std::arithmetic auto x) : m_x(static_cast<num_t>(x)) {}

  template<std::arithmetic T>
  constexpr explicit operator T () const {
      return static_cast<T>(m_x);
  }

  constexpr std::partial_ordering operator<=>(const unreal_t& oth) const {
      if (definitelyLesserThan(*this, oth)) {
      return std::partial_ordering::less;
    }

    if(definitelyGreaterThan(*this, oth)) {
      return std::partial_ordering::greater;
    }

    if(approximatelyEqual(*this, oth)) {
      return std::partial_ordering::equivalent;
    }

    return std::partial_ordering::unordered;
  }

  constexpr bool operator==(const unreal_t& oth) const = default;
  constexpr bool operator!=(const unreal_t& oth) const = default;
  constexpr bool operator<=(const unreal_t& oth) const = default;
  constexpr bool operator>=(const unreal_t& oth) const = default;
  constexpr bool operator< (const unreal_t& oth) const = default;
  constexpr bool operator> (const unreal_t& oth) const = default;

  static constexpr bool
  approximatelyEqual(unreal_t lhs, unreal_t rhs) {
    return std::abs(lhs.m_x - rhs.m_x) <=
           ((std::abs(lhs.m_x) < std::abs(rhs.m_x) ? std::abs(rhs.m_x)
                                                   : std::abs(lhs.m_x)) *
            Epsilon);
  }

  static constexpr bool
  essentialyEqual(unreal_t lhs, unreal_t rhs) {
    return std::abs(lhs.m_x - rhs.m_x) <=
           ((std::abs(lhs.m_x) > std::abs(rhs.m_x) ? std::abs(rhs.m_x)
                                                   : std::abs(lhs.m_x)) *
            Epsilon);
  }

  static constexpr bool
  definitelyGreaterThan(unreal_t lhs, unreal_t rhs) {
    if(lhs.m_x == std::numeric_limits<num_t>::infinity() || rhs.m_x == -std::numeric_limits<num_t>::infinity()) {
      return true;
    }
    return (lhs.m_x - rhs.m_x) >
           ((std::abs(lhs.m_x) < std::abs(rhs.m_x) ? std::abs(rhs.m_x)
                                                   : std::abs(lhs.m_x)) *
            Epsilon);
  }

  static constexpr bool
  definitelyLesserThan(unreal_t lhs, unreal_t rhs) {
    if(lhs.m_x == -std::numeric_limits<num_t>::infinity() || rhs.m_x == std::numeric_limits<num_t>::infinity()) {
      return true;
    }
    return (rhs.m_x - lhs.m_x) >
           ((std::abs(lhs.m_x) < std::abs(rhs.m_x) ? std::abs(rhs.m_x)
                                                   : std::abs(lhs.m_x)) *
            Epsilon);
  }

  constexpr unreal_t& operator+= (unreal_t oth) {
    m_x += oth.m_x;
    return *this;
  }

  constexpr unreal_t& operator-= (unreal_t oth) {
    m_x -= oth.m_x;
    return *this;
  }

  constexpr unreal_t& operator*= (unreal_t oth) {
    m_x *= oth.m_x;
    return *this;
  }

  constexpr unreal_t& operator/= (unreal_t oth) {
    m_x /= oth.m_x;
    return *this;
  }

  constexpr unreal_t operator+ (unreal_t oth) const {
    return unreal_t{m_x} += oth.m_x;
  }

  constexpr unreal_t operator- (unreal_t oth) const {
    return unreal_t{m_x} -= oth.m_x;
  }

  constexpr unreal_t operator* (unreal_t oth) const {
    return unreal_t{m_x} *= oth.m_x;
  }

  constexpr unreal_t operator/ (unreal_t oth) const {
    return unreal_t{m_x} /= oth.m_x;
  }

  constexpr unreal_t operator- () const {
    return unreal_t{-m_x};
  }

  constexpr num_t getVal() const {return m_x;}
};

namespace std {

constexpr unreal_t abs(unreal_t x) {
  return unreal_t{abs(x.getVal())};
}

constexpr unreal_t sqrt(unreal_t x) {
  return unreal_t{sqrt(x.getVal())};
}

}

template<typename OStream>
OStream& operator<<(OStream& out, const unreal_t x) {
  return out << x.getVal();
}

template<typename IStream>
IStream& operator>>(IStream& in, unreal_t& x) {
  unreal_t::num_t xv;
  in >> xv;
  x = xv;
  return in;
}

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUG__)
#pragma GCC diagnostic pop
#else
#endif

#endif // REAL_HPP
