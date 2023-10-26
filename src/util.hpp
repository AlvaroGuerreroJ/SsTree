#pragma once

#include "params.h"
#include "Point.h"

struct Pair {
    Point point;
    NType distance;

    Pair(const Point& p, NType d) : point(p), distance(d) {}
};

struct Comparator {
    bool operator()(const Pair& a, const Pair& b) const {
        return a.distance < b.distance; // max-heap basado en distancia
    }
};

template <typename It>
double variance(It begin, It end) {
    size_t count = 0;
    double sum = 0;

    for (auto it = begin; it != end; it++) {
        sum += *it;
        count++;
    }

    double avg = sum / count;

    double variance = 0;
    for (auto it = begin; it != end; it++) {
        double d = *it - avg;
        variance += d * d;
    }
    variance /= count;

    return variance;
}
