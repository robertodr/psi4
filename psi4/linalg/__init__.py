#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Linear algebra utilities"""

from ..core.linalg import *
from .factory import Matrix_, Tensor3, Tensor4, Vector_
from .utils import (make_random_tensor_2d, make_tensor_2d_from_block, name_unary_test)

# NOTE this is rather awkward, but I haven't found any other way to get the
# help for all the stuff defined in core.linalg otherwise...
# yapf: disable
__all__ = [
      "Vector_" , "Vector_F", "Vector_D", "Vector_CD"
    , "Matrix_", "Matrix_F", "Matrix_D", "Matrix_CD"
    , "Tensor3", "Tensor3_F", "Tensor3_D", "Tensor3_CD"
    , "Tensor4", "Tensor4_F", "Tensor4_D", "Tensor4_CD"
    , "Operation"
    , "gemm"
    , "doublet"
    , "gemv"
    , "zeros_like"
    , "ones_like"
    , "full_like"
    , "real"
    , "imag"
    , "conj"
    , "cholesky"
    , "eig"
    , "make_random_tensor_2d"
    , "make_tensor_2d_from_block"
    , "name_unary_test"
    , "Vector3"
    , "dot", "norm", "cross", "normalize", "distance", "angle", "perp_unit", "rotate"
]
# yapf: enable