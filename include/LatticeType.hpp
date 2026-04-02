#pragma once
#include <vector>
#include <cmath>

template <typename T>
struct lattice {
    size_t Nx, Ny, Nz;
    size_t ySkip, zSkip;
    std::vector<T> field;

    lattice(size_t const x, size_t const y, size_t const z): Nx(x), Ny(y), Nz(z), ySkip(x), zSkip(x * y), field(x*y*z) {}

    const T &operator()(size_t const x, size_t const y, size_t const z) const {
        return field[x + ySkip * y + zSkip * z];
    }

    T &operator()(size_t const x, size_t const y, size_t const z) {
        return field[x + ySkip * y + zSkip * z];
    }

    const T &operator[](size_t const x) const {
        return field[x];
    }

    T &operator[](size_t const x) {
        return field[x];
    }

    void zero() {
        for (auto &element: field) {
            element = 0;
        }
    }

    T mean() const {
        T total{};
        #pragma omp parallel for
        for (int i = 0; i < field.size(); ++i) {
            total += field[i];
        }
        return total / field.size();
    }

    lattice<T> &operator+=(const lattice<T> &rhs) {
        #pragma omp parallel for
        for (size_t i = 0; i < field.size(); ++i)
            field[i] += rhs.field[i];
        return *this;
    }

    lattice<T> &operator-=(const lattice<T> &rhs) {
        #pragma omp parallel for
        for (size_t i = 0; i < field.size(); ++i)
            field[i] -= rhs.field[i];
        return *this;
    }

    lattice<T> &operator*=(const lattice<T> &rhs) {
        #pragma omp parallel for
        for (size_t i = 0; i < field.size(); ++i)
            field[i] *= rhs.field[i];
        return *this;
    }

    lattice<T> &operator/=(const lattice<T> &rhs) {
        #pragma omp parallel for
        for (size_t i = 0; i < field.size(); ++i)
            field[i] /= rhs.field[i];
        return *this;
    }



    lattice<T> &operator+=(const T &rhs) {
        #pragma omp parallel for
        for (size_t i = 0; i < field.size(); ++i)
            field[i] += rhs;
        return *this;
    }

    lattice<T> &operator-=(const T &rhs) {
        #pragma omp parallel for
        for (size_t i = 0; i < field.size(); ++i)
            field[i] -= rhs;
        return *this;
    }

    lattice<T> &operator*=(const T &rhs) {
        #pragma omp parallel for
        for (size_t i = 0; i < field.size(); ++i)
            field[i] *= rhs;
        return *this;
    }

    lattice<T> &operator/=(const T &rhs) {
        #pragma omp parallel for
        for (size_t i = 0; i < field.size(); ++i)
            field[i] /= rhs;
        return *this;
    }


};

template <typename T>
lattice<T> operator+(lattice<T> lhs, const lattice<T> &rhs) {
    lhs += rhs;
    return lhs;
}

template <typename T>
lattice<T> operator-(lattice<T> lhs, const lattice<T> &rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator*(lattice<T> lhs, const lattice<T> &rhs) {
    lhs *= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator/(lattice<T> lhs, const lattice<T> &rhs) {
    lhs /= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator+(lattice<T> lhs, const T &rhs) {
    lhs += rhs;
    return lhs;
}

template <typename T>
lattice<T> operator-(lattice<T> lhs, const T &rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator*(lattice<T> lhs, const T &rhs) {
    lhs *= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator/(lattice<T> lhs, const T &rhs) {
    lhs /= rhs;
    return lhs;
}