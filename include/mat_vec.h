/**
 * Copyright (C) 2023-2024 Simon Hafner - All Rights Reserved  
 * Institute of Flight System Dynamics -- Technische Universitaet Muenchen
 * Contact: simon.hafner@tum.de  
 */

#ifndef MATRIX_H
#define MATRIX_H

// #define MAXROWS 6
/* Specific to project and must be adapted*/

#define EPS_F 1.19e-07f

#define ABS(X) (((X) >= (float)0.0) ? (X) : -(X)) /* Absolutbetrag von X */

#define MIN(X, Y) (((X) <= (Y)) ? (X) : (Y)) /* Minimum of X and Y */
#define MAX(X, Y) (((X) >= (Y)) ? (X) : (Y)) /* Maximum of X and Y  */

#define M_P_ENTRY(mat, i, j) ((mat)->matrix[((j) * (mat)->rows) + (i)])
#define M_ENTRY(mat, i, j) ((mat).matrix[((j) * (mat).rows) + (i)])

// #include <math.h>

/*Variable definitions*/
/*(C) S. Myschik, OCT 2005*/

/**
 * @brief Matrix struct containing matrix, num rows, num columns
 *
 */
typedef struct
{
    unsigned short rows;
    unsigned short columns;
    float *matrix;
} matrix_t;

/**
 * @brief Constant Matrix struct containing matrix, num rows, num columns
 *
 */
typedef struct
{
    const unsigned short rows;
    const unsigned short columns;
    const float *matrix;
} matrix_c_t;

/*(C) S. Hafner. OCT 2023*/

/**
 * @brief Vector struct containing matrix, num rows
 *
 */
typedef struct
{
    unsigned short rows;
    float *vector;
} vector_t;

/**
 * @brief Constant Vector struct containing matrix, num rows
 *
 */
typedef struct
{
    const unsigned short rows;
    const float *vector;
} vector_c_t;

/**
 * @brief Vector struct containing matrix, num rows
 *
 */
typedef struct
{
    unsigned short rows;
    unsigned short *vector;
} vector_u_t;

/**
 * @brief Vector struct containing matrix, num rows
 *
 */
typedef struct
{
    unsigned short rows;
    short *vector;
} vector_i_t;

/*Function definitions*/
/*(C) S. Myschik, OCT 2005*/

/*Matrix operating  functions*/
#ifdef TEST_CONFIG
void mat_create(matrix_t *matr, unsigned short rows, unsigned short cols);
#endif

void eye(matrix_t *matr, unsigned short elements);

// void mat_print(matrix_t matr);

// void mat_clear(matrix_t *matr);

void mat_zero(matrix_t *const matr);

void mat_copy(const matrix_t *const matr, matrix_t *const res_matr);

unsigned int mat_mult(const matrix_t *const matr1, const matrix_t *const matr2, matrix_t *const res_matr);

void mat_mult_scalar(matrix_t *const matr, const float scalar);

/*(C) S. Hafner. OCT 2023*/
#ifdef TEST_CONFIG
void vec_create(vector_t *vec, unsigned short rows);
void vec_create_u(vector_u_t *vec, unsigned short rows);
void print_matrix(const matrix_t *const matr);
#endif

void vec_zero(vector_t *const vec);

float vec_max(const vector_t *const vec);

void vec_scalar(vector_t *const vec, const float scalar);

unsigned int vec_pivot(const vector_t *const vec, const vector_u_t *const pivot, vector_t *res_vec);

unsigned int vec_pivot_c(const vector_c_t *const vec, const vector_u_t *const pivot, vector_t *res_vec);

unsigned int vec_pivot_rev(const vector_t *const vec, const vector_u_t *const pivot, vector_t *res_vec);

unsigned int vec_inv(vector_t *const vec);

unsigned int vec_T_vec(const vector_t *const vec1, const vector_t *const vec2, float *res);

unsigned int vec_sum(const vector_t *const vec1, const vector_t *const vec2, vector_t *res_vec);

unsigned int vec_diff(const vector_t *const vec1, const vector_t *const vec2, vector_t *res_vec);

void vec_raw_copy(const float *vec, vector_t *const res_vec);

void vec_copy(const vector_t *const vec, vector_t *const res_vec);

// TODO vec_clear
void mat_raw_copy(const float *matr, matrix_t *const res_matr);

void mat_diag(const matrix_t *const matr, vector_t *res_vec);

float mat_max_diag(const matrix_t *const matr);

unsigned short mat_max_diag_pivot(const matrix_t *const matr, const vector_u_t *const p, unsigned short start);

unsigned int mat_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec);

unsigned int mat_t_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec);

unsigned int mat_t_vec_c(const matrix_c_t *const matr, const vector_t *const vec, vector_t *res_vec);

unsigned int diag_vec(const vector_t *const diag, const vector_t *const vec, vector_t *res_vec);

unsigned int diag_vec_c(const vector_c_t *const diag, const vector_t *const vec, vector_t *res_vec);

unsigned int diag_inv_vec(const vector_t *const diag, const vector_t *const vec, vector_t *res_vec);

unsigned int lower_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec);

unsigned int lower_t_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec);

unsigned int lower_t_vec_c(const matrix_t *const matr, const vector_c_t *const vec, vector_t *res_vec);

unsigned int lower_inv(const matrix_t *const matr, matrix_t *res_matr);

unsigned int lower_inv_LDLT(const matrix_t *const matr, matrix_t *res_matr);

float det_two(const matrix_t *const matr, const unsigned short row_1, const unsigned short row_2, const unsigned short col_1, const unsigned short col_2);

float det_three(const matrix_t *const matr, const unsigned short row_1, const unsigned short row_2, const unsigned short row_3, const unsigned short col_1, const unsigned short col_2, const unsigned short col_3);

unsigned int mat_gram(const matrix_t *const matr, matrix_t *res_matr);

unsigned int mat_gram_weighted(const matrix_t *const matr, const vector_t *const W, matrix_t *res_matr);

unsigned int mat_gram_weighted_c(const matrix_c_t *const matr, const vector_c_t *const W, matrix_t *res_matr);

unsigned int mat_gram_T_weighted(const matrix_t *const matr, const vector_t *const W, matrix_t *res_matr);

unsigned int mat_mult_diag_c(matrix_t *const matr, const vector_c_t *const W);

unsigned int lower_T_lower(const matrix_t *const matr, matrix_t *res_matr);

unsigned int solve_lower_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec);

unsigned int solve_lower_T_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec);

unsigned int solve_LDL_lower_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec);

unsigned int solve_LDL_lower_vec_c(const matrix_t *const matr, const vector_c_t *const vec, vector_t *res_vec);

unsigned int solve_LDL_lower_T_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec);

#endif /* MATRIX */