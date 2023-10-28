#ifndef ENGINE_PHYSCONSTANTS_HPP
#define ENGINE_PHYSCONSTANTS_HPP

#include "units.hpp"

namespace phys {
namespace consts {

constexpr const auto G = 6.67430e-11_H * 1._m * 1._m / 1._kg / 1._kg;
constexpr const auto k = 1.38064e-23_J / 1_K;
constexpr const auto Dalton = 1.660e-27_kg;
} // namespace consts

namespace config {

constexpr const num_t DT_NORMALIZER = 4e-5;
constexpr const num_t FLUCTATION_DEGREE = 1e3;

} // namespace config
} // namespace phys

#endif /* ENGINE_PHYSCONSTANTS_HPP */
