/**
 * Copyright (C) 2023-2024 Simon Hafner - All Rights Reserved  
 * Institute of Flight System Dynamics -- Technische Universitaet Muenchen
 * Contact: simon.hafner@tum.de  
 */

#ifndef SYM_MAT_INV_H
#define SYM_MAT_INV_H

#include "mat_vec.h"

unsigned int sym_inv_two(const matrix_t *const matr, matrix_t *res_matr);

unsigned int sym_inv_three(const matrix_t *const matr, matrix_t *res_matr);

unsigned int sym_inv_four(const matrix_t *const matr, matrix_t *res_matr);

unsigned int sym_inv_five(const matrix_t *const matr, matrix_t *res_matr);

// unsigned int sym_inv_six(matrix_t matr, matrix_t *res_matr);

#endif /*SYM_MAT_INV_H*/