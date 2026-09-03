#pragma once
#include <cstddef>
#include <vector>
#include <array>

template <typename T>
struct lattice {
    std::size_t Nx, Ny, Nz;
    std::vector<T> field;

    lattice(): Nx(0), Ny(0), Nz(0), field(0) {};
    lattice(std::size_t const x, std::size_t const y, std::size_t const z): Nx(x), Ny(y), Nz(z), field(x*y*z) {}

    T* data() {
        return field.data();
    }

    T const* data() const {
        return field.data();
    }

    T const &operator()(std::size_t const x, std::size_t const y, std::size_t const z) const {
        return field[x + Nx * (y + Ny * z)];
    }

    T &operator()(std::size_t const x, std::size_t const y, std::size_t const z) {
        return field[x + Nx * ( y + Ny * z)];
    }

    T const &operator[](std::size_t const x) const {
        return field[x];
    }

    T &operator[](std::size_t const x) {
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
        for (std::size_t i = 0; i < field.size(); ++i) {
            total += field[i];
        }
        return total / field.size();
    }

    T size() const {
        return field.size();
    }

    lattice<T> &operator+=(lattice<T> const &rhs) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] += rhs.field[i];
        return *this;
    }

    lattice<T> &operator-=(lattice<T> const &rhs) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] -= rhs.field[i];
        return *this;
    }

    lattice<T> &operator*=(lattice<T> const &rhs) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] *= rhs.field[i];
        return *this;
    }

    lattice<T> &operator/=(lattice<T> const &rhs) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] /= rhs.field[i];
        return *this;
    }


    lattice<T> &operator+=(T const &rhs) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] += rhs;
        return *this;
    }

    lattice<T> &operator-=(T const &rhs) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] -= rhs;
        return *this;
    }

    lattice<T> &operator*=(T const &rhs) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] *= rhs;
        return *this;
    }

    lattice<T> &operator/=(T const &rhs) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] /= rhs;
        return *this;
    }


};

template <typename T>
lattice<T> operator+(lattice<T> lhs, lattice<T> const &rhs) {
    lhs += rhs;
    return lhs;
}

template <typename T>
lattice<T> operator-(lattice<T> lhs, lattice<T> const &rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator*(lattice<T> lhs, lattice<T> const &rhs) {
    lhs *= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator/(lattice<T> lhs, lattice<T> const &rhs) {
    lhs /= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator+(lattice<T> lhs, T const &rhs) {
    lhs += rhs;
    return lhs;
}

template <typename T>
lattice<T> operator-(lattice<T> lhs, T const &rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator*(lattice<T> lhs, T const &rhs) {
    lhs *= rhs;
    return lhs;
}

template <typename T>
lattice<T> operator/(lattice<T> lhs, T const &rhs) {
    lhs /= rhs;
    return lhs;
}