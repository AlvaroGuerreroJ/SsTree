#pragma once

#include <cmath>
#include <vector>

#include "Point.h"

auto vector_norm(std::vector<double> const& v) -> double;
auto vector_unit(std::vector<double> const& v) -> std::vector<double>;
auto vector_scale(std::vector<double> const& v, double c) -> std::vector<double>;
auto vector_add(std::vector<double> const& v1, std::vector<double> const& v2)
    -> std::vector<double>;
auto vector_sub(std::vector<double> const& v1, std::vector<double> const& v2)
    -> std::vector<double>;

auto point_to_vd(Point const& p) -> std::vector<double>;
