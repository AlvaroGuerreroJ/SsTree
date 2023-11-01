#include <cassert>
#include <cmath>
#include <vector>

#include "pvector.hpp"

auto vector_norm(std::vector<double> const& v) -> double
{
    double r = 0;
    for (double d : v)
    {
        r += d * d;
    }
    return std::sqrt(r);
}

auto vector_unit(std::vector<double> const& v) -> std::vector<double>
{
    double n = vector_norm(v);
    std::vector<double> ret;
    for (double d : v)
    {
        ret.push_back(d / n);
    }

    return ret;
}

auto vector_scale(std::vector<double> const& v, double c) -> std::vector<double>
{
    std::vector<double> ret;
    for (double d : v)
    {
        ret.push_back(c * d);
    }

    return ret;
}

std::vector<double> vector_add(std::vector<double> const& v1, std::vector<double> const& v2)
{
    assert(v1.size() == v2.size());

    std::vector<double> ret(v1.size());

    for (size_t i = 0; i < v1.size(); i++)
    {
        ret[i] = v1[i] + v2[i];
    }

    return ret;
}

/* Substracts v2 from v1.
 *
 * i.e. [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2], …]
 */
auto vector_sub(std::vector<double> const& v1, std::vector<double> const& v2)
    -> std::vector<double>
{
    assert(v1.size() == v2.size());

    std::vector<double> ret(v1.size());

    for (size_t i = 0; i < v1.size(); i++)
    {
        ret[i] = v1[i] - v2[i];
    }

    return ret;
}

auto point_to_vd(Point const& p) -> std::vector<double>
{
    return {p.begin(), p.end()};
}
