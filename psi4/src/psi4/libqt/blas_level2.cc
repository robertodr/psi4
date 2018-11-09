/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

/**!
 * \file
 * \brief Interface to all BLAS routines
 * \ingroup QT
 *
 * Autogenerated by Rob Parrish on 1/24/2011
 *
 */

#include "blas_level2.h"

#include <stdexcept>

#ifdef USING_LAPACK_MKL
#include "mkl_cblas.h"
#else
#include "cblas.h"
#endif

namespace psi {
void C_DGBMV(char trans, int m, int n, int kl, int ku, double alpha, double* a, int lda, double* x, int incx,
             double beta, double* y, int incy) {
    if (m == 0 || n == 0) return;
    if (trans == 'N' || trans == 'n')
        trans = 'T';
    else if (trans == 'T' || trans == 't')
        trans = 'N';
    else
        throw std::invalid_argument("C_DGBMV trans argument is invalid.");
    ::F_DGBMV(&trans, &n, &m, &ku, &kl, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void C_DGEMV(char trans, int m, int n, double alpha, double* a, int lda, double* x, int incx, double beta, double* y,
             int incy) {
    if (m == 0 || n == 0) return;
    if (trans == 'N' || trans == 'n')
        trans = 'T';
    else if (trans == 'T' || trans == 't')
        trans = 'N';
    else
        throw std::invalid_argument("C_DGEMV trans argument is invalid.");
    ::F_DGEMV(&trans, &n, &m, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void C_DGER(int m, int n, double alpha, double* x, int incx, double* y, int incy, double* a, int lda) {
    if (m == 0 || n == 0) return;
    ::F_DGER(&n, &m, &alpha, y, &incy, x, &incx, a, &lda);
}

void C_DSBMV(char uplo, int n, int k, double alpha, double* a, int lda, double* x, int incx, double beta, double* y,
             int incy) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DSBMV uplo argument is invalid.");
    ::F_DSBMV(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void C_DSPMV(char uplo, int n, double alpha, double* ap, double* x, int incx, double beta, double* y, int incy) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DSPMV uplo argument is invalid.");
    ::F_DSPMV(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
}

void C_DSPR(char uplo, int n, double alpha, double* x, int incx, double* ap) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DSPR uplo argument is invalid.");
    ::F_DSPR(&uplo, &n, &alpha, x, &incx, ap);
}

void C_DSPR2(char uplo, int n, double alpha, double* x, int incx, double* y, int incy, double* ap) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DSPR2 uplo argument is invalid.");
    ::F_DSPR2(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
}

void C_DSYMV(char uplo, int n, double alpha, double* a, int lda, double* x, int incx, double beta, double* y,
             int incy) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DSYMV uplo argument is invalid.");
    ::F_DSYMV(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void C_DSYR(char uplo, int n, double alpha, double* x, int incx, double* a, int lda) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DSYR uplo argument is invalid.");
    ::F_DSYR(&uplo, &n, &alpha, x, &incx, a, &lda);
}

void C_DSYR2(char uplo, int n, double alpha, double* x, int incx, double* y, int incy, double* a, int lda) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DSYR2 uplo argument is invalid.");
    ::F_DSYR2(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

void C_DTBMV(char uplo, char trans, char diag, int n, int k, double* a, int lda, double* x, int incx) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DTBMV uplo argument is invalid.");
    if (trans == 'N' || trans == 'n')
        trans = 'T';
    else if (trans == 'T' || trans == 't')
        trans = 'N';
    else
        throw std::invalid_argument("C_DTBMV trans argument is invalid.");
    ::F_DTBMV(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

void C_DTBSV(char uplo, char trans, char diag, int n, int k, double* a, int lda, double* x, int incx) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DTBSV uplo argument is invalid.");
    if (trans == 'N' || trans == 'n')
        trans = 'T';
    else if (trans == 'T' || trans == 't')
        trans = 'N';
    else
        throw std::invalid_argument("C_DTBSV trans argument is invalid.");
    ::F_DTBSV(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

void C_DTPMV(char uplo, char trans, char diag, int n, double* ap, double* x, int incx) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DTPMV uplo argument is invalid.");
    if (trans == 'N' || trans == 'n')
        trans = 'T';
    else if (trans == 'T' || trans == 't')
        trans = 'N';
    else
        throw std::invalid_argument("C_DTPMV trans argument is invalid.");
    ::F_DTPMV(&uplo, &trans, &diag, &n, ap, x, &incx);
}

void C_DTPSV(char uplo, char trans, char diag, int n, double* ap, double* x, int incx) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DTPSV uplo argument is invalid.");
    if (trans == 'N' || trans == 'n')
        trans = 'T';
    else if (trans == 'T' || trans == 't')
        trans = 'N';
    else
        throw std::invalid_argument("C_DTPSV trans argument is invalid.");
    ::F_DTPSV(&uplo, &trans, &diag, &n, ap, x, &incx);
}

void C_DTRMV(char uplo, char trans, char diag, int n, double* a, int lda, double* x, int incx) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DTRMV uplo argument is invalid.");
    if (trans == 'N' || trans == 'n')
        trans = 'T';
    else if (trans == 'T' || trans == 't')
        trans = 'N';
    else
        throw std::invalid_argument("C_DTRMV trans argument is invalid.");
    ::F_DTRMV(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

void C_DTRSM(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b,
             int ldb) {
    if (m == 0 || n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DTRSM uplo argument is invalid.");
    if (side == 'L' || side == 'L')
        side = 'R';
    else if (side == 'R' || side == 'r')
        side = 'L';
    else
        throw std::invalid_argument("C_DTRSM side argument is invalid.");
    ::F_DTRSM(&side, &uplo, &transa, &diag, &n, &m, &alpha, a, &lda, b, &ldb);
}

void C_DTRSV(char uplo, char trans, char diag, int n, double* a, int lda, double* x, int incx) {
    if (n == 0) return;
    if (uplo == 'U' || uplo == 'u')
        uplo = 'L';
    else if (uplo == 'L' || uplo == 'l')
        uplo = 'U';
    else
        throw std::invalid_argument("C_DTRSV uplo argument is invalid.");
    if (trans == 'N' || trans == 'n')
        trans = 'T';
    else if (trans == 'T' || trans == 't')
        trans = 'N';
    else
        throw std::invalid_argument("C_DTRSV trans argument is invalid.");
    ::F_DTRSV(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}
}  // namespace psi
