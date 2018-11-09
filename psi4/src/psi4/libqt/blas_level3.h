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

/*!
 * \file
 * \brief Function signatures for wrappers to BLAS level 3 subroutines
 * \ingroup QT
 */

#pragma once

#include "psi4/pragma.h"

namespace psi {
/**
 *  Purpose
 *  =======
 *
 *  DGEMM  performs one of the matrix-matrix operations
 *
 *     C := alpha*op( A )*op( B ) + beta*C,
 *
 *  where  op( X ) is one of
 *
 *     op( X ) = X   or   op( X ) = X',
 *
 *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
 *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
 *
 *  Arguments
 *  ==========
 *
 *  TRANSA - CHARACTER*1.
 *           On entry, TRANSA specifies the form of op( A ) to be used in
 *           the matrix multiplication as follows:
 *
 *              TRANSA = 'N' or 'n',  op( A ) = A.
 *
 *              TRANSA = 'T' or 't',  op( A ) = A'.
 *
 *              TRANSA = 'C' or 'c',  op( A ) = A'.
 *
 *           Unchanged on exit.
 *
 *  TRANSB - CHARACTER*1.
 *           On entry, TRANSB specifies the form of op( B ) to be used in
 *           the matrix multiplication as follows:
 *
 *              TRANSB = 'N' or 'n',  op( B ) = B.
 *
 *              TRANSB = 'T' or 't',  op( B ) = B'.
 *
 *              TRANSB = 'C' or 'c',  op( B ) = B'.
 *
 *           Unchanged on exit.
 *
 *  M      - INTEGER.
 *           On entry,  M  specifies  the number  of rows  of the  matrix
 *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry,  N  specifies the number  of columns of the matrix
 *           op( B ) and the number of columns of the matrix C. N must be
 *           at least zero.
 *           Unchanged on exit.
 *
 *  K      - INTEGER.
 *           On entry,  K  specifies  the number of columns of the matrix
 *           op( A ) and the number of rows of the matrix op( B ). K must
 *           be at least  zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
 *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
 *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
 *           part of the array  A  must contain the matrix  A,  otherwise
 *           the leading  k by m  part of the array  A  must contain  the
 *           matrix A.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
 *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
 *           least  max( 1, k ).
 *           Unchanged on exit.
 *
 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
 *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
 *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
 *           part of the array  B  must contain the matrix  B,  otherwise
 *           the leading  n by k  part of the array  B  must contain  the
 *           matrix B.
 *           Unchanged on exit.
 *
 *  LDB    - INTEGER.
 *           On entry, LDB specifies the first dimension of B as declared
 *           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
 *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
 *           least  max( 1, n ).
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
 *           supplied as zero then C need not be set on input.
 *           Unchanged on exit.
 *
 *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
 *           Before entry, the leading  m by n  part of the array  C must
 *           contain the matrix  C,  except when  beta  is zero, in which
 *           case C need not be set on entry.
 *           On exit, the array  C  is overwritten by the  m by n  matrix
 *           ( alpha*op( A )*op( B ) + beta*C ).
 *
 *  LDC    - INTEGER.
 *           On entry, LDC specifies the first dimension of C as declared
 *           in  the  calling  (sub)  program.   LDC  must  be  at  least
 *           max( 1, m ).
 *           Unchanged on exit.
 *
 *
 *  Level 3 Blas routine.
 *
 *  -- Written on 8-February-1989.
 *     Jack Dongarra, Argonne National Laboratory.
 *     Iain Duff, AERE Harwell.
 *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     Sven Hammarling, Numerical Algorithms Group Ltd.
 *
 *
 *     .. External Functions ..
 */
PSI_API
void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb,
             double beta, double* c, int ldc);

/**
 *  Purpose
 *  =======
 *
 *  DSYMM  performs one of the matrix-matrix operations
 *
 *     C := alpha*A*B + beta*C,
 *
 *  or
 *
 *     C := alpha*B*A + beta*C,
 *
 *  where alpha and beta are scalars,  A is a symmetric matrix and  B and
 *  C are  m by n matrices.
 *
 *  Arguments
 *  ==========
 *
 *  SIDE   - CHARACTER*1.
 *           On entry,  SIDE  specifies whether  the  symmetric matrix  A
 *           appears on the  left or right  in the  operation as follows:
 *
 *              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
 *
 *              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
 *
 *           Unchanged on exit.
 *
 *  UPLO   - CHARACTER*1.
 *           On  entry,   UPLO  specifies  whether  the  upper  or  lower
 *           triangular  part  of  the  symmetric  matrix   A  is  to  be
 *           referenced as follows:
 *
 *              UPLO = 'U' or 'u'   Only the upper triangular part of the
 *                                  symmetric matrix is to be referenced.
 *
 *              UPLO = 'L' or 'l'   Only the lower triangular part of the
 *                                  symmetric matrix is to be referenced.
 *
 *           Unchanged on exit.
 *
 *  M      - INTEGER.
 *           On entry,  M  specifies the number of rows of the matrix  C.
 *           M  must be at least zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the number of columns of the matrix C.
 *           N  must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
 *           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
 *           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
 *           the array  A  must contain the  symmetric matrix,  such that
 *           when  UPLO = 'U' or 'u', the leading m by m upper triangular
 *           part of the array  A  must contain the upper triangular part
 *           of the  symmetric matrix and the  strictly  lower triangular
 *           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
 *           the leading  m by m  lower triangular part  of the  array  A
 *           must  contain  the  lower triangular part  of the  symmetric
 *           matrix and the  strictly upper triangular part of  A  is not
 *           referenced.
 *           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
 *           the array  A  must contain the  symmetric matrix,  such that
 *           when  UPLO = 'U' or 'u', the leading n by n upper triangular
 *           part of the array  A  must contain the upper triangular part
 *           of the  symmetric matrix and the  strictly  lower triangular
 *           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
 *           the leading  n by n  lower triangular part  of the  array  A
 *           must  contain  the  lower triangular part  of the  symmetric
 *           matrix and the  strictly upper triangular part of  A  is not
 *           referenced.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
 *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
 *           least  max( 1, n ).
 *           Unchanged on exit.
 *
 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
 *           Before entry, the leading  m by n part of the array  B  must
 *           contain the matrix B.
 *           Unchanged on exit.
 *
 *  LDB    - INTEGER.
 *           On entry, LDB specifies the first dimension of B as declared
 *           in  the  calling  (sub)  program.   LDB  must  be  at  least
 *           max( 1, m ).
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
 *           supplied as zero then C need not be set on input.
 *           Unchanged on exit.
 *
 *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
 *           Before entry, the leading  m by n  part of the array  C must
 *           contain the matrix  C,  except when  beta  is zero, in which
 *           case C need not be set on entry.
 *           On exit, the array  C  is overwritten by the  m by n updated
 *           matrix.
 *
 *  LDC    - INTEGER.
 *           On entry, LDC specifies the first dimension of C as declared
 *           in  the  calling  (sub)  program.   LDC  must  be  at  least
 *           max( 1, m ).
 *           Unchanged on exit.
 *
 *
 *  Level 3 Blas routine.
 *
 *  -- Written on 8-February-1989.
 *     Jack Dongarra, Argonne National Laboratory.
 *     Iain Duff, AERE Harwell.
 *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     Sven Hammarling, Numerical Algorithms Group Ltd.
 *
 *
 *     .. External Functions ..
 */
PSI_API
void C_DSYMM(char side, char uplo, int m, int n, double alpha, double* a, int lda, double* b, int ldb, double beta,
             double* c, int ldc);

/**
 *  Purpose
 *  =======
 *
 *  DTRMM  performs one of the matrix-matrix operations
 *
 *     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
 *
 *  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
 *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
 *
 *     op( A ) = A   or   op( A ) = A'.
 *
 *  Arguments
 *  ==========
 *
 *  SIDE   - CHARACTER*1.
 *           On entry,  SIDE specifies whether  op( A ) multiplies B from
 *           the left or right as follows:
 *
 *              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
 *
 *              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
 *
 *           Unchanged on exit.
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the matrix A is an upper or
 *           lower triangular matrix as follows:
 *
 *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 *
 *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 *
 *           Unchanged on exit.
 *
 *  TRANSA - CHARACTER*1.
 *           On entry, TRANSA specifies the form of op( A ) to be used in
 *           the matrix multiplication as follows:
 *
 *              TRANSA = 'N' or 'n'   op( A ) = A.
 *
 *              TRANSA = 'T' or 't'   op( A ) = A'.
 *
 *              TRANSA = 'C' or 'c'   op( A ) = A'.
 *
 *           Unchanged on exit.
 *
 *  DIAG   - CHARACTER*1.
 *           On entry, DIAG specifies whether or not A is unit triangular
 *           as follows:
 *
 *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *
 *              DIAG = 'N' or 'n'   A is not assumed to be unit
 *                                  triangular.
 *
 *           Unchanged on exit.
 *
 *  M      - INTEGER.
 *           On entry, M specifies the number of rows of B. M must be at
 *           least zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the number of columns of B.  N must be
 *           at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
 *           zero then  A is not referenced and  B need not be set before
 *           entry.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
 *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
 *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
 *           upper triangular part of the array  A must contain the upper
 *           triangular matrix  and the strictly lower triangular part of
 *           A is not referenced.
 *           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
 *           lower triangular part of the array  A must contain the lower
 *           triangular matrix  and the strictly upper triangular part of
 *           A is not referenced.
 *           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
 *           A  are not referenced either,  but are assumed to be  unity.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
 *           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
 *           then LDA must be at least max( 1, n ).
 *           Unchanged on exit.
 *
 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
 *           Before entry,  the leading  m by n part of the array  B must
 *           contain the matrix  B,  and  on exit  is overwritten  by the
 *           transformed matrix.
 *
 *  LDB    - INTEGER.
 *           On entry, LDB specifies the first dimension of B as declared
 *           in  the  calling  (sub)  program.   LDB  must  be  at  least
 *           max( 1, m ).
 *           Unchanged on exit.
 *
 *
 *  Level 3 Blas routine.
 *
 *  -- Written on 8-February-1989.
 *     Jack Dongarra, Argonne National Laboratory.
 *     Iain Duff, AERE Harwell.
 *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     Sven Hammarling, Numerical Algorithms Group Ltd.
 *
 *
 *     .. External Functions ..
 */
PSI_API
void C_DTRMM(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b,
             int ldb);

/**
 *  Purpose
 *  =======
 *
 *  DSYRK  performs one of the symmetric rank k operations
 *
 *     C := alpha*A*A' + beta*C,
 *
 *  or
 *
 *     C := alpha*A'*A + beta*C,
 *
 *  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
 *  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
 *  in the second case.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On  entry,   UPLO  specifies  whether  the  upper  or  lower
 *           triangular  part  of the  array  C  is to be  referenced  as
 *           follows:
 *
 *              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
 *                                  is to be referenced.
 *
 *              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
 *                                  is to be referenced.
 *
 *           Unchanged on exit.
 *
 *  TRANS  - CHARACTER*1.
 *           On entry,  TRANS  specifies the operation to be performed as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
 *
 *              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
 *
 *              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry,  N specifies the order of the matrix C.  N must be
 *           at least zero.
 *           Unchanged on exit.
 *
 *  K      - INTEGER.
 *           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
 *           of  columns   of  the   matrix   A,   and  on   entry   with
 *           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
 *           of rows of the matrix  A.  K must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
 *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
 *           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
 *           part of the array  A  must contain the matrix  A,  otherwise
 *           the leading  k by n  part of the array  A  must contain  the
 *           matrix A.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
 *           then  LDA must be at least  max( 1, n ), otherwise  LDA must
 *           be at least  max( 1, k ).
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry, BETA specifies the scalar beta.
 *           Unchanged on exit.
 *
 *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
 *           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
 *           upper triangular part of the array C must contain the upper
 *           triangular part  of the  symmetric matrix  and the strictly
 *           lower triangular part of C is not referenced.  On exit, the
 *           upper triangular part of the array  C is overwritten by the
 *           upper triangular part of the updated matrix.
 *           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
 *           lower triangular part of the array C must contain the lower
 *           triangular part  of the  symmetric matrix  and the strictly
 *           upper triangular part of C is not referenced.  On exit, the
 *           lower triangular part of the array  C is overwritten by the
 *           lower triangular part of the updated matrix.
 *
 *  LDC    - INTEGER.
 *           On entry, LDC specifies the first dimension of C as declared
 *           in  the  calling  (sub)  program.   LDC  must  be  at  least
 *           max( 1, n ).
 *           Unchanged on exit.
 *
 *
 *  Level 3 Blas routine.
 *
 *  -- Written on 8-February-1989.
 *     Jack Dongarra, Argonne National Laboratory.
 *     Iain Duff, AERE Harwell.
 *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     Sven Hammarling, Numerical Algorithms Group Ltd.
 *
 *
 *     .. External Functions ..
 */
PSI_API
void C_DSYRK(char uplo, char trans, int n, int k, double alpha, double* a, int lda, double beta, double* c, int ldc);

/**
 *  Purpose
 *  =======
 *
 *  DSYR2K  performs one of the symmetric rank 2k operations
 *
 *     C := alpha*A*B' + alpha*B*A' + beta*C,
 *
 *  or
 *
 *     C := alpha*A'*B + alpha*B'*A + beta*C,
 *
 *  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
 *  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
 *  matrices in the second case.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On  entry,   UPLO  specifies  whether  the  upper  or  lower
 *           triangular  part  of the  array  C  is to be  referenced  as
 *           follows:
 *
 *              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
 *                                  is to be referenced.
 *
 *              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
 *                                  is to be referenced.
 *
 *           Unchanged on exit.
 *
 *  TRANS  - CHARACTER*1.
 *           On entry,  TRANS  specifies the operation to be performed as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
 *                                        beta*C.
 *
 *              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
 *                                        beta*C.
 *
 *              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
 *                                        beta*C.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry,  N specifies the order of the matrix C.  N must be
 *           at least zero.
 *           Unchanged on exit.
 *
 *  K      - INTEGER.
 *           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
 *           of  columns  of the  matrices  A and B,  and on  entry  with
 *           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
 *           of rows of the matrices  A and B.  K must be at least  zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
 *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
 *           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
 *           part of the array  A  must contain the matrix  A,  otherwise
 *           the leading  k by n  part of the array  A  must contain  the
 *           matrix A.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
 *           then  LDA must be at least  max( 1, n ), otherwise  LDA must
 *           be at least  max( 1, k ).
 *           Unchanged on exit.
 *
 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
 *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
 *           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
 *           part of the array  B  must contain the matrix  B,  otherwise
 *           the leading  k by n  part of the array  B  must contain  the
 *           matrix B.
 *           Unchanged on exit.
 *
 *  LDB    - INTEGER.
 *           On entry, LDB specifies the first dimension of B as declared
 *           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
 *           then  LDB must be at least  max( 1, n ), otherwise  LDB must
 *           be at least  max( 1, k ).
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry, BETA specifies the scalar beta.
 *           Unchanged on exit.
 *
 *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
 *           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
 *           upper triangular part of the array C must contain the upper
 *           triangular part  of the  symmetric matrix  and the strictly
 *           lower triangular part of C is not referenced.  On exit, the
 *           upper triangular part of the array  C is overwritten by the
 *           upper triangular part of the updated matrix.
 *           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
 *           lower triangular part of the array C must contain the lower
 *           triangular part  of the  symmetric matrix  and the strictly
 *           upper triangular part of C is not referenced.  On exit, the
 *           lower triangular part of the array  C is overwritten by the
 *           lower triangular part of the updated matrix.
 *
 *  LDC    - INTEGER.
 *           On entry, LDC specifies the first dimension of C as declared
 *           in  the  calling  (sub)  program.   LDC  must  be  at  least
 *           max( 1, n ).
 *           Unchanged on exit.
 *
 *
 *  Level 3 Blas routine.
 *
 *
 *  -- Written on 8-February-1989.
 *     Jack Dongarra, Argonne National Laboratory.
 *     Iain Duff, AERE Harwell.
 *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     Sven Hammarling, Numerical Algorithms Group Ltd.
 *
 *
 *     .. External Functions ..
 */
PSI_API
void C_DSYR2K(char uplo, char trans, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta,
              double* c, int ldc);
}  // namespace psi
