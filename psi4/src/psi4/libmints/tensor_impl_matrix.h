/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#pragma once

#include <memory>
#include <string>

#include <xtensor/xadapt.hpp>

#include "dimension.h"
#include "matrix.h"
#include "tensor_impl.h"

namespace psi {
template <typename T, size_t Rank, typename Enable>
class Tensor;

template <typename T>
using Matrix_ = Tensor<T, 2, std::enable_if_t<detail::is_tensorisable_v<T>>>;

template <typename T>
using SharedMatrix_ = std::shared_ptr<Matrix_<T>>;

namespace detail {
template <typename T>
struct RankDependentImpl<Matrix_<T>> {
    void set_symmetric(size_t h, size_t i, size_t j, T val) {
        static_cast<Matrix_<T>*>(this)->store_[h](i, j) = static_cast<Matrix_<T>*>(this)->store_[h](j, i) = val;
    }
    void set_symmetric(size_t i, size_t j, T val) { this->set_symmetric(0, i, j, val); }

    /*! C++ string representation of type */
    std::string cxxClassName() const noexcept { return "Matrix<" + detail::Type2String<T>::full() + ">"; }
    /*! Python string representation of type */
    static std::string pyClassName() noexcept { return "Matrix_" + detail::Type2String<T>::suffix(); }

    /// Returns the dimension array for the rows of a rank-2 tensor aka a matrix
    const Dimension& rowspi() const { return static_cast<const Matrix_<T>*>(this)->axes_dimpi_[0]; }
    /*! Returns the number of rows in given irrep
     *  \param[in] h
     */
    size_t rows(size_t h = 0) const { return static_cast<const Matrix_<T>*>(this)->axes_dimpi_[0][h]; }

    /// Returns the dimension array for the rows of a rank-2 tensor aka a matrix
    const Dimension& colspi() const { return static_cast<const Matrix_<T>*>(this)->axes_dimpi_[1]; }
    /*! Returns the number of columns in given irrep
     *  \param[in] h
     */
    size_t cols(size_t h = 0) const { return static_cast<const Matrix_<T>*>(this)->axes_dimpi_[1][h]; }
};
}  // namespace detail

/*! Conversion from Matrix to Matrix_<double>
 *  \param[in] m
 */
template <typename T = double>
auto transmute(const Matrix& m) -> Matrix_<T> {
    auto mat = Matrix_<T>(m.name(), m.rowspi(), m.colspi(), m.symmetry());
    // Adapt memory owned by matrix class
    for (int h = 0; h < m.nirrep(); ++h) {
        std::vector<std::size_t> shape = {m.rowdim(h), m.coldim(h)};
        mat[h] = xt::adapt(m.get_pointer(h), m.rowdim(h) * m.coldim(h), xt::no_ownership(), shape);
    }
    return mat;
}

/*! Conversion from SharedMatrix to SharedMatrix_<double>
 *  \param[in] m
 */
template <typename T = double>
auto transmute(const std::shared_ptr<Matrix>& m) -> SharedMatrix_<T> {
    auto mat = std::make_shared<Matrix_<T>>(m->name(), m->rowspi(), m->colspi(), m->symmetry());
    // Adapt memory owned by matrix class
    for (int h = 0; h < m->nirrep(); ++h) {
        std::vector<std::size_t> shape = {m->rowdim(h), m->coldim(h)};
        mat->block(h) = xt::adapt(m->get_pointer(h), m->rowdim(h) * m->coldim(h), xt::no_ownership(), shape);
    }
    return mat;
}

/*! Conversion from Matrix_<double> to Matrix
 *  \param[in] m
 */
template <typename T = double>
auto transmute(const Matrix_<T>& m) -> Matrix {
    auto mat = Matrix(m.label(), m.rowspi(), m.colspi(), m.symmetry());
    // Element-by-element copy
    for (int h = 0; h < m.nirrep(); ++h) {
        for (auto i = 0; i < m.rows(h); ++i) {
            for (auto j = 0; i < m.cols(h); ++j) {
                mat.set(h, i, j, m.get(h, i, j));
            }
        }
    }
    return mat;
}

/*! Conversion from SharedMatrix_<double> to SharedMatrix
 *  \param[in] v
 */
template <typename T = double>
auto transmute(const SharedMatrix_<T>& m) -> SharedMatrix {
    auto mat = std::make_shared<Matrix>(m->label(), m->rowspi(), m->colspi(), m->symmetry());
    // Element-by-element copy
    for (int h = 0; h < m->nirrep(); ++h) {
        for (auto i = 0; i < m->rows(h); ++i) {
            for (auto j = 0; i < m->cols(h); ++j) {
                mat->set(h, i, j, m->get(h, i, j));
            }
        }
    }
    return mat;
}
}  // namespace psi