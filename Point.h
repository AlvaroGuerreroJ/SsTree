#ifndef POINT_H
#define POINT_H

#include "params.h"
#include <vector>
#include <cmath>


class Point {
private:
    std::vector<NType> coordinates;

public:
    Point() {}
    Point(size_t size) : coordinates(size) {}
    Point(std::initializer_list<NType> init) : coordinates(init) {}

    // Acceso a las coordenadas
    const NType& operator[](size_t index) const {
        return coordinates[index];
    }
    NType& operator[](size_t index) {
        return coordinates[index];
    }

    // Operadores de comparación
    bool operator==(const Point& other) const {
        for (size_t i = 0; i < coordinates.size(); ++i) {
            if (coordinates[i] != other.coordinates[i]) {
                return false;
            }
        }
        return true;
    }
    bool operator!=(const Point& other) const {
        return !(*this == other);
    }

    // Operaciones
    Point operator+(const Point& other) const {
        Point result(dim());
        for (size_t i = 0; i < dim(); ++i) {
            result[i] = coordinates[i] + other.coordinates[i];
        }
        return result;
    }
    Point operator-(const Point& other) const {
        Point result(dim());
        for (size_t i = 0; i < dim(); ++i) {
            result[i] = coordinates[i] - other.coordinates[i];
        }
        return result;
    }
    Point operator*(NType scalar) const {
        Point result(dim());
        for (size_t i = 0; i < dim(); ++i) {
            result[i] = coordinates[i] * scalar;
        }
        return result;
    }
    Point operator/(NType scalar) const {
        Point result(dim());
        for (size_t i = 0; i < dim(); ++i) {
            result[i] = coordinates[i] / scalar;
        }
        return result;
    }

    // Operadores de asignación
    Point& operator+=(const Point& other) {
        for (size_t i = 0; i < dim(); ++i) {
            coordinates[i] += other.coordinates[i];
        }
        return *this;
    }
    Point& operator-=(const Point& other) {
        for (size_t i = 0; i < dim(); ++i) {
            coordinates[i] -= other.coordinates[i];
        }
        return *this;
    }
    Point& operator*=(NType scalar) {
        for (size_t i = 0; i < dim(); ++i) {
            coordinates[i] *= scalar;
        }
        return *this;
    }
    Point& operator/=(NType scalar) {
        for (size_t i = 0; i < dim(); ++i) {
            coordinates[i] /= scalar;
        }
        return *this;
    }

    // Funciones auxiliares
    size_t dim() const {
        return coordinates.size();
    }
    NType norm() const {
        NType result = 0;
        for (size_t i = 0; i < dim(); ++i) {
            result += coordinates[i] * coordinates[i];
        }
        return sqrt(result);
    }

    // Iteradores para recorrer las coordenadas
    using iterator = typename std::vector<NType>::iterator;
    using const_iterator = typename std::vector<NType>::const_iterator;
    iterator begin() {
        return coordinates.begin();
    }
    const_iterator begin() const {
        return coordinates.begin();
    }
    iterator end() {
        return coordinates.end();
    }
    const_iterator end() const {
        return coordinates.end();
    }

    // Funciones de utilidad
    friend std::ostream& operator<<(std::ostream& os, const Point& p){
        os << "(";
        for (auto it = p.begin(); it != p.end(); ++it) {
            os << *it;
            if (it + 1 != p.end()) {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
};

NType distance(const Point& a, const Point& b) {
    if (a.dim() != b.dim()) {
        throw std::runtime_error("Los puntos deben tener la misma dimensión");
    }
    NType sum = 0;
    for (auto it1 = a.begin(), it2 = b.begin(); it1 != a.end(); ++it1, ++it2) {
        sum = sum + (*it1 - *it2) * (*it1 - *it2);
    }
    return NType::sqrt(sum);
}
NType manhattanDistance(const Point& a, const Point& b) {
    if (a.dim() != b.dim()) {
        throw std::runtime_error("Los puntos deben tener la misma dimensión");
    }
    NType sum = 0;
    for (auto it1 = a.begin(), it2 = b.begin(); it1 != a.end(); ++it1, ++it2) {
        sum = sum + NType::abs(*it1 - *it2);
    }
    return sum;
}
NType chebyshevDistance(const Point& a, const Point& b) {
    if (a.dim() != b.dim()) {
        throw std::runtime_error("Los puntos deben tener la misma dimensión");
    }
    NType maxDiff = 0;
    for (auto it1 = a.begin(), it2 = b.begin(); it1 != a.end(); ++it1, ++it2) {
        maxDiff = NType::max(maxDiff, NType::abs(*it1 - *it2));
    }
    return maxDiff;
}
NType minkowskiDistance(const Point& a, const Point& b, int p) {
    if (a.dim() != b.dim()) {
        throw std::runtime_error("Los puntos deben tener la misma dimensión");
    }
    if (p <= 0) {
        throw std::runtime_error("El valor de p debe ser positivo");
    }
    NType sum = 0;
    for (auto it1 = a.begin(), it2 = b.begin(); it1 != a.end(); ++it1, ++it2) {
        sum = sum + NType::pow(NType::abs(*it1 - *it2), p);
    }
    return NType::pow(sum, 1.0 / p);
}


#endif // !POINT_H