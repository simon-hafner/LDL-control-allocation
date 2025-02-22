/**
 * Copyright (C) 2023-2024 Simon Hafner - All Rights Reserved  
 * Institute of Flight System Dynamics -- Technische Universitaet Muenchen
 * Contact: simon.hafner@tum.de  
 */

#ifndef CONTROL_ALLOCATION_LDL_H
#define CONTROL_ALLOCATION_LDL_H

#define TOL 1e-5f

#include "mat_vec.h"

unsigned int CA_full_LDLT_six(matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, unsigned short *const rank, vector_t *u_alloc, const unsigned short rank_max);

unsigned int sim_CA_full_LDLT_six(const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max);

unsigned int sim_CA_full_LDLT_four(const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max);

unsigned int sim_CA_piv_LDLT_six(const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max);

unsigned int sim_CA_piv_LDLT_four(const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max);

unsigned int allocation_LDLT_six(const matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, unsigned short *const rank, vector_t *u_alloc);

unsigned int allocation_LDLT_four(const matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, unsigned short *const rank, vector_t *u_alloc);

unsigned int allocation_piv_LDLT_six(const matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, vector_u_t *p, unsigned short *const rank, vector_t *u_alloc);

unsigned int allocation_piv_LDLT_four(const matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, vector_u_t *p, unsigned short *const rank, vector_t *u_alloc);

unsigned int CA_down_LDLT_six(const unsigned short saturated, matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, unsigned short *const rank, vector_t *u_alloc, unsigned short rank_max);

unsigned int sim_CA_down_LDLT_six(const unsigned short saturated, const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max);

unsigned int sim_CA_down_LDLT_four(const unsigned short saturated, const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max);

unsigned int full_rank_Cholesky(const matrix_t *const A, matrix_t *L, unsigned short *const rank);

unsigned int full_rank_LDLT(const matrix_t *const A, matrix_t *L, vector_t *d, unsigned short *const rank, const float tol, const unsigned short rank_max);

unsigned int pivoted_LDLT(const matrix_t *const A, matrix_t *L, vector_t *d, vector_u_t *p, unsigned short *const rank, const float tol, const unsigned short rank_max);

unsigned int rank_one_LDLT_down(vector_t v, matrix_t *L, vector_t *d, unsigned short *const rank, const float tol, const unsigned short rank_max);

#endif /* CONTROL_ALLOCATION_LDL_H */
