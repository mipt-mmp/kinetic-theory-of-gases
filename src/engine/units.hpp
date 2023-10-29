#ifndef ENGINE_UNITS_HPP
#define ENGINE_UNITS_HPP
#include "geometry.hpp"
#include "real.hpp"
#include <compare>
#include <concepts>
#include <cstdint>

#define PHYS_UNITS_PROVIDE_LITERALS // FIXME: Compile option

namespace phys {

const std::size_t UniverseDim = 3;

namespace detail {
class UnitBase {}; // Needed for fast abstraction
} // namespace detail

template <class T>
concept SomeUnit = std::derived_from<T, detail::UnitBase>;

template <typename Num, int Meter = 0, int Sec = 0, int KG = 0, int K = 0, int Mole = 0,
          int Ampere = 0, int Candela = 0>
class Unit : public detail::UnitBase {
    Num num_;

public:
    constexpr Unit()
        : num_(0) {}

    constexpr explicit Unit(const Num& t)
        : num_(t) {}

    constexpr Num inner() {
        return num_;
    }

    constexpr Num& operator*() {
        return num_;
    }

    constexpr const Num& operator*() const {
        return num_;
    }

    constexpr Num* operator->() {
        return &num_;
    }

    constexpr const Num* operator->() const {
        return &num_;
    }

    constexpr Unit& operator+=(const Unit& rhs) {
        num_ += rhs.num_;
        return *this;
    }

    constexpr Unit operator+(const Unit& rhs) const {
        return Unit{*this} += rhs;
    }

    constexpr Unit& operator-=(const Unit& rhs) {
        num_ -= rhs.num_;
        return *this;
    }

    constexpr Unit operator-(const Unit& rhs) const {
        return Unit{*this} -= rhs;
    }

    constexpr Unit operator-() const {
        return Unit{-num_};
    }

    constexpr Unit& operator*=(const Num& num) {
        num_ *= num;
        return *this;
    }

    constexpr Unit& operator/=(const Num& num) {
        num_ /= num;
        return *this;
    }

    constexpr Unit operator/(const Num& num) {
        return Unit{*this} /= num;
    }

    constexpr std::partial_ordering operator<=>(const Unit& rhs) const {
        return num_ <=> rhs.num_;
    }

    constexpr bool operator==(const Unit& rhs) const = default;
    constexpr bool operator!=(const Unit& rhs) const = default;
    constexpr bool operator<=(const Unit& rhs) const = default;
    constexpr bool operator>=(const Unit& rhs) const = default;
    constexpr bool operator<(const Unit& rhs) const = default;
    constexpr bool operator>(const Unit& rhs) const = default;
};

template <typename Num, int Meter1, int Sec1, int KG1, int K1, int Mole1, int Ampere1, int Candela1,
          int Meter2, int Sec2, int KG2, int K2, int Mole2, int Ampere2, int Candela2>
constexpr auto operator*(const Unit<Num, Meter1, Sec1, KG1, K1, Mole1, Ampere1, Candela1>& lhs,
                         const Unit<Num, Meter2, Sec2, KG2, K2, Mole2, Ampere2, Candela2>& rhs)
    -> Unit<Num, Meter1 + Meter2, Sec1 + Sec2, KG1 + KG2, K1 + K2, Mole1 + Mole2, Ampere1 + Ampere2,
            Candela1 + Candela2> {
    return Unit<Num, Meter1 + Meter2, Sec1 + Sec2, KG1 + KG2, K1 + K2, Mole1 + Mole2,
                Ampere1 + Ampere2, Candela1 + Candela2>{*lhs * *rhs};
}

template <typename Num, int Meter, int Sec, int KG, int K, int Mole, int Ampere, int Candela>
constexpr auto operator*(const Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela>& lhs,
                         const Num& rhs) -> Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela> {
    return Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela>{*lhs} *= rhs;
}

template <typename Num, int Meter, int Sec, int KG, int K, int Mole, int Ampere, int Candela>
constexpr auto operator*(const Num& lhs,
                         const Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela>& rhs)
    -> Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela> {
    return Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela>{*rhs} *= lhs;
}

template <typename Num, int Meter1, int Sec1, int KG1, int K1, int Mole1, int Ampere1, int Candela1,
          int Meter2, int Sec2, int KG2, int K2, int Mole2, int Ampere2, int Candela2>
constexpr auto operator/(const Unit<Num, Meter1, Sec1, KG1, K1, Mole1, Ampere1, Candela1>& lhs,
                         const Unit<Num, Meter2, Sec2, KG2, K2, Mole2, Ampere2, Candela2>& rhs)
    -> Unit<Num, Meter1 - Meter2, Sec1 - Sec2, KG1 - KG2, K1 - K2, Mole1 - Mole2, Ampere1 - Ampere2,
            Candela1 - Candela2> {
    return Unit<Num, Meter1 - Meter2, Sec1 - Sec2, KG1 - KG2, K1 - K2, Mole1 - Mole2,
                Ampere1 - Ampere2, Candela1 - Candela2>{*lhs / *rhs};
}

template <typename Num, int Meter, int Sec, int KG, int K, int Mole, int Ampere, int Candela>
constexpr auto sqrt(const Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela>& u)
    -> Unit<Num, Meter / 2, Sec / 2, KG / 2, K / 2, Mole / 2, Ampere / 2, Candela / 2> {
    static_assert(Meter % 2 == 0, "Bad sqrt");
    static_assert(Sec % 2 == 0, "Bad sqrt");
    static_assert(KG % 2 == 0, "Bad sqrt");
    static_assert(K % 2 == 0, "Bad sqrt");
    static_assert(Mole % 2 == 0, "Bad sqrt");
    static_assert(Ampere % 2 == 0, "Bad sqrt");
    static_assert(Candela % 2 == 0, "Bad sqrt");
    return Unit<Num, Meter / 2, Sec / 2, KG / 2, K / 2, Mole / 2, Ampere / 2, Candela / 2>{
        std::sqrt(*u)};
}

template <SomeUnit T>
struct UnitTraits {
    static auto sqrt(const T& t) -> decltype(phys::sqrt(t)) {
        return phys::sqrt(t);
    }
};

using num_t = unreal_t;

using LengthVal = Unit<num_t, 1>; // m

using SquareVal = Unit<num_t, 2>; // m^2

using VolumeVal = Unit<num_t, UniverseDim>; // m^3

using VelocityVal = Unit<num_t, 1, -1>; // m * s^-1

using AccelerationVal = Unit<num_t, 1, -2>; // m * s^-2

using ForceVal = Unit<num_t, 1, -2, 1>; // kg * m * s^-2

using PressureVal = decltype(ForceVal{} * LengthVal{} / VolumeVal{}); // kg * m^-1 * s^-2

using TimeVal = Unit<num_t, 0, 1, 0>; // s

using MassVal = Unit<num_t, 0, 0, 1>; // kg

using EnergyVal = Unit<num_t, 2, -2, 1>; // kg * m ^ 2 * s^-2

using MassMomentVal = Unit<num_t, 1, 0, 1>; // kg * m

using ImpulseVal = Unit<num_t, 1, -1, 1>; // kg * m * s^-1

using ImpulseMomentVal = Unit<num_t, 2, -1, 1>; // kg * m^2 * s^-1

using TemperatureVal = Unit<num_t, 0, 0, 0, 1>; // kg * m^2 * s^-1

using Length = LengthVal;
using Square = SquareVal;
using Volume = VolumeVal;
using Pressure = PressureVal;
using Temperature = TemperatureVal;

template <SomeUnit T>
using Vector = geom::Vector<T, UniverseDim, UnitTraits>;

using Position = Vector<LengthVal>;
using Distance = Vector<LengthVal>;
using MassMoment = Vector<MassMomentVal>;

using Velocity = Vector<VelocityVal>;
using Acceleration = Vector<AccelerationVal>;
using Force = Vector<ForceVal>;
using Impulse = Vector<ImpulseVal>;
using ImpulseMoment = Vector<ImpulseMomentVal>;

using Time = TimeVal;
using Mass = MassVal;
using Energy = EnergyVal;

constexpr const num_t Zepto = 1e-21;
constexpr const num_t Atto = 1e-18;
constexpr const num_t Femto = 1e-15;
constexpr const num_t Pico = 1e-12;
constexpr const num_t Nano = 1e-9;
constexpr const num_t Micro = 1e-6;
constexpr const num_t Milli = 1e-3;

constexpr const num_t Kilo = 1e3;
constexpr const num_t Mega = 1e6;
constexpr const num_t Giga = 1e9;
constexpr const num_t Tera = 1e12;
constexpr const num_t Peta = 1e15;
constexpr const num_t Exa = 1e18;
constexpr const num_t Zetta = 1e21;
constexpr const num_t Yotta = 1e24;
constexpr const num_t Ronna = 1e27;
constexpr const num_t Quetta = 1e30;

constexpr static const std::pair<num_t, const char*> suffixes[] = {
    {Zepto, "z"}, {Atto, "a"}, {Femto, "f"}, {Pico, "p"},  {Nano, "n"},  {Micro, "u"},
    {Milli, "m"}, {1e0, ""},   {Kilo, "K"},  {Mega, "M"},  {Giga, "G"},  {Tera, "T"},
    {Peta, "P"},  {Exa, "E"},  {Zetta, "Z"}, {Yotta, "Y"}, {Ronna, "R"}, {Quetta, "Q"},
};

template <SomeUnit T>
std::pair<num_t, const char*> getSuff(const T& unit) {
    if (std::abs(*unit) / suffixes[0].first < 1.) {
        return suffixes[0];
    }
    for (size_t i = 1; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
        if (std::abs(*unit) / suffixes[i].first < 1.) {
            return suffixes[i - 1];
        }
    }
    return suffixes[sizeof(suffixes) / sizeof(suffixes[0]) - 1];
}

template <SomeUnit T, SomeUnit U>
auto operator*(const Vector<T>& lhs, const U& rhs) -> Vector<decltype(lhs[0] * rhs)> {
    Vector<decltype(lhs[0] * rhs)> ans{};
    for (size_t i = 0; i < UniverseDim; ++i) {
        ans[i] = lhs[i] * rhs;
    }
    return ans;
}

template <SomeUnit T, SomeUnit U>
auto MulByElement(const Vector<T>& lhs, const Vector<U>& rhs) -> Vector<decltype(lhs[0] * rhs[0])> {
    Vector<decltype(lhs[0] * rhs[0])> ans{};
    for (size_t i = 0; i < UniverseDim; ++i) {
        ans[i] = lhs[i] * rhs[i];
    }
    return ans;
}

template <SomeUnit T, typename U>
auto operator/(const Vector<T>& lhs, const U& rhs) -> Vector<decltype(lhs[0] / rhs)> {
    Vector<decltype(lhs[0] / rhs)> ans{};
    for (size_t i = 0; i < UniverseDim; ++i) {
        ans[i] = lhs[i] / rhs;
    }
    return ans;
}

template <SomeUnit T>
auto Normalize(const Vector<T>& t) {
    assert(*t.Len() > 0.);
    return t / t.Len();
}

template <SomeUnit T>
T Trace(const Vector<T>& t) {
    T ans{};
    for (size_t i = 0; i < UniverseDim; ++i) {
        ans += t[i];
    }
    return ans;
}

template <SomeUnit T>
Vector<T> randomSphere() {
    Vector<T> v;
    for (size_t i = 0; i < UniverseDim; ++i) {
        v[i] = T{num_t{rand()}};
    }
    v /= *v.Len();
    return v;
}

template <SomeUnit T>
Vector<T> randomInCube(Vector<T> mx) {
    Vector<T> v;
    for (size_t i = 0; i < UniverseDim; ++i) {
        v[i] = mx[i] * num_t{rand()};
    }
    v /= num_t{RAND_MAX};
    return v;
}

static inline num_t randomShift() {
    return num_t{(rand() - RAND_MAX / 2)} / num_t{RAND_MAX};
}

template <SomeUnit T, SomeUnit U>
Unit<num_t> GetSinusBetween(const Vector<T>& lhs, const Vector<U>& rhs) {
    auto lens_product = lhs.Len() * rhs.Len();
    auto cos_alpha = (lhs, rhs) / lens_product;
    auto sin_alpha = Unit<num_t>{1} - cos_alpha * cos_alpha;

    return sin_alpha;
}

#define PRINT_UNIT(Unit, Name)                \
    do {                                      \
        if constexpr (Unit != 0) {            \
            name += Name;                     \
            if constexpr (Unit != 1) {        \
                name += "^";                  \
                name += std::to_string(Unit); \
            }                                 \
        }                                     \
    } while (0)

template <typename T>
struct UnitFormatter;

template <typename Num, int Meter, int Sec, int KG, int K, int Mole, int Ampere, int Candela>
struct UnitFormatter<Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela>> {
    static constexpr std::string str() {
        std::string name = "";
        PRINT_UNIT(Meter, "m");
        PRINT_UNIT(Sec, "s");
        PRINT_UNIT(KG, "kg");
        PRINT_UNIT(K, "K");
        PRINT_UNIT(Mole, "mol");
        PRINT_UNIT(Ampere, "A");
        PRINT_UNIT(Candela, "cd");
        return name;
    }
};

#undef PRINT_UNIT

template <>
struct UnitFormatter<ForceVal> {
    static constexpr std::string str() {
        return "H";
    }
};

template <>
struct UnitFormatter<EnergyVal> {
    static constexpr std::string str() {
        return "J";
    }
};

template <>
struct UnitFormatter<PressureVal> {
    static constexpr std::string str() {
        return "Pa";
    }
};

} // namespace phys

#if defined(PHYS_UNITS_PROVIDE_LITERALS)
#    define PHYS_UNIT_LITERAL(U, Literal)                              \
        constexpr phys::U operator""_##Literal(long double x) {        \
            return phys::U{x};                                         \
        }                                                              \
        constexpr phys::U operator""_##Literal(unsigned long long x) { \
            return phys::U{x};                                         \
        }

PHYS_UNIT_LITERAL(Length, m)
PHYS_UNIT_LITERAL(Time, sec)
PHYS_UNIT_LITERAL(Mass, kg)
PHYS_UNIT_LITERAL(ForceVal, H)
PHYS_UNIT_LITERAL(EnergyVal, J)
PHYS_UNIT_LITERAL(Temperature, K)
#endif

template <class IStream, typename Num, int Meter, int Sec, int KG, int K, int Mole, int Ampere,
          int Candela>
IStream& operator>>(IStream& in, phys::Unit<Num, Meter, Sec, KG, K, Mole, Ampere, Candela>& u) {
    return in >> *u;
}

template <class OStream, phys::SomeUnit T>
OStream& operator<<(OStream& out, const T& u) {

    auto [div, pref] = phys::getSuff(u);
    out << *u / div << ' ' << pref << phys::UnitFormatter<T>::str().c_str();

    return out;
}

#endif /* ENGINE_UNITS_HPP */
