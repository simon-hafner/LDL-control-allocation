/**
 * Copyright (C) 2023-2024 Simon Hafner - All Rights Reserved  
 * Institute of Flight System Dynamics -- Technische Universitaet Muenchen
 * Contact: simon.hafner@tum.de  
 */

#ifdef TEST_CONFIG
#include <stdlib.h>
#include <stdio.h>
#endif
#include <math.h>
#include "mat_vec.h"

#ifdef TEST_CONFIG
/**
 * @brief Creates dynamiocally a new zero n x m matrix (only for debugging purposes)
 *
 * @param matr Pointer to the new matrix
 * @param rows Number of rows n
 * @param cols Number of columns n
 */
void mat_create(matrix_t *matr, unsigned short rows, unsigned short cols)
{

    matr->columns = cols; /*  Save number of columns */
    matr->rows = rows;    /* Save number of rows  */

    matr->matrix = (float *)malloc((float)cols * (float)rows * sizeof(float)); /*  Allocate memory for each row */
    for (unsigned short n = 0U; n < rows; n++)                                 /*  Fill all elements with 0 */
    {
        for (unsigned short m = 0U; m < cols; m++)
        {
            (M_P_ENTRY(matr, n, m)) = 0.0f;
        }
    }
}
#endif

/**
 * @brief Sets a n x n matrix to a unity/identity matrix
 * @copyright Myschik, OCT-2005
 *
 * @param matr Pointer to the identity matrix
 * @param elements Number of elements n
 */
void eye(matrix_t *matr, unsigned short elements)
{
    // TODO Check for square
    for (unsigned short n = 0U; n < elements; n++) /* Fill all elements with 0 */
    {
        for (unsigned short m = 0U; m < elements; m++)
        {
            if (m == n)
            {
                (M_P_ENTRY(matr, n, n)) = 1.0f;
            }
            else
            {
                (M_P_ENTRY(matr, n, m)) = 0.0f;
            }
        }
    }
}

/**
 * @brief Prints content of m x n Matrix  (for debugging purpose)
 * @copyright Myschik, OCT-2005
 *
 * @param matr Matrix to be printed
 */
// void mat_print(matrix_t matr)
// {
//     int no_of_cols = matr.columns;
//     int no_of_rows = matr.rows;
//     printf("\n");

//     for (int n = 0;n < no_of_rows; n++)
//     {
//         for (int m = 0; m < no_of_cols; m++)
//         {
//             printf("%.10f \t",matr.matrix[n + (rows * m)]);
//         }
//         printf("\n");
//     }
// }

/**
 * @brief Fills every element of matrix with 0
 * @copyright Myschik, OCT-2005
 *
 * @param matr Pointer to the matrix to be overwritten
 */
void mat_zero(matrix_t *const matr)
{
    unsigned short no_of_rows = matr->rows; /*  get number of rows  */
    unsigned short no_of_cols = matr->columns;

    for (unsigned short n = 0U; n < no_of_rows; n++)
    {
        for (unsigned short m = 0U; m < no_of_cols; m++)
        {
            M_P_ENTRY(matr, n, m) = 0.0f;
        }
    }
}

/**
 * @brief Elementwise copy of mxn matrices
 * @copyright Myschik, OCT-2005
 *
 * @param matr Matrix to be copied (original matrix)
 * @param res_matr Pointer to target matrix
 */
void mat_copy(const matrix_t *const matr, matrix_t *const res_matr)
{
    unsigned short no_of_cols = matr->columns; /* Maximum Number of Columns */
    unsigned short no_of_rows = matr->rows;    /* Maximum Number of Rows */

    // TODO Check for right dimensions
    for (unsigned short n = 0U; n < no_of_rows; n++)
    {
        for (unsigned short m = 0U; m < no_of_cols; m++)
        {
            M_P_ENTRY(res_matr, n, m) = M_P_ENTRY(matr, n, m); /* Print each Column element in one row*/
        }
    }
}

/**
 * @brief Multiplicates two m x n matrices  A*B = C
 * @copyright Myschik, OCT-2005
 *
 * @param matr1 Matrix A (n x m)
 * @param matr2 Matrix B (m x p)
 * @param res_matr Resulting Matrix C (n x p)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_mult(const matrix_t *const matr1, const matrix_t *const matr2, matrix_t *const res_matr)
{
    unsigned int status = 0;
    unsigned short left_col_size;  /* Number of Columns of Matrix A (left) */
    unsigned short left_row_size;  /* Number of Rows of Matrix A (left) */
    unsigned short right_col_size; /* Number of Rows of Matrix B (right) */
    unsigned short right_row_size; /* Number of Rows of Matrix B (right) */

    left_col_size = matr1->columns; /* get size of left matrix */
    left_row_size = matr1->rows;

    right_col_size = matr2->columns; /* get size of right matrix */
    right_row_size = matr2->rows;

    if (right_row_size != left_col_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        mat_zero(res_matr);
        for (unsigned short h = 0U; h < left_row_size; h++) /* left row counter               */
        {
            for (unsigned short j = 0U; j < right_col_size; j++) /*  right column counter          */
            {
                for (unsigned short i = 0U; i < left_col_size; i++) /* left col counter               */
                {
                    M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, h, j) + /* calculate result     */
                                                (M_P_ENTRY(matr1, h, i) * M_P_ENTRY(matr2, i, j));
                }
            }
        }
        /* return status 0 (OK)         */
    }
    return (status);
}

/**
 * @brief Elementwise multiplication of scalar with matrix $B = \sigma * A$
 *
 * @param matr Input matrix A and output B (A is overwritten)
 * @param scalar Scalar sigma
 */
void mat_mult_scalar(matrix_t *const matr, const float scalar)
{
    unsigned short no_of_cols = matr->columns; /* Maximum Number of Columns */
    unsigned short no_of_rows = matr->rows;    /* Maximum Number of Rows */

    for (unsigned short n = 0U; n < no_of_rows; n++)
    {
        for (unsigned short m = 0U; m < no_of_cols; m++)
        {
            M_P_ENTRY(matr, n, m) = M_P_ENTRY(matr, n, m) * scalar; /* Multiply each element with scalar */
        }
    }
}

/***************************************************************************/
/*(C) S. Hafner. OCT 2023*/
/***************************************************************************/
#ifdef TEST_CONFIG
/**
 * @brief vec_create creates a new zero vector struct with a specific number of rows
 * (only for debugging puposes)
 *
 * @param vec pointer to the vector struct
 * @param rows number of rows in the new vector
 */
void vec_create(vector_t *vec, unsigned short rows)
{
    vec->rows = rows; /* Save number of rows  */

    vec->vector = (float *)malloc((float)rows * sizeof(float)); /*  Allocate memory for each row */

    for (unsigned short n = 0U; n < rows; n++) /*  Fill all elements with 0 */
    {
        (vec->vector[n]) = 0.0f;
    }
}

void vec_create_u(vector_u_t *vec, unsigned short rows)
{
    vec->rows = rows; /* Save number of rows  */

    vec->vector = (unsigned short *)malloc((float)rows * sizeof(unsigned short)); /*  Allocate memory for each row */

    for (unsigned short n = 0U; n < rows; n++) /*  Fill all elements with 0 */
    {
        (vec->vector[n]) = 0U;
    }
}

void print_matrix(const matrix_t *const matr)
{
    for (unsigned short n = 0U; n < matr->rows; n++)
    {
        for (unsigned short m = 0U; m < matr->columns; m++)
        {
            printf("%.8f ", M_P_ENTRY(matr, n, m));
        }
        printf("\n");
    }
}

#endif

/**
 * @brief Set all elements of a vector vec to zero
 *
 * @param vec Vector to be set to zero
 */
void vec_zero(vector_t *const vec)
{
    for (unsigned short i = 0U; i < vec->rows; i++)
    {
        vec->vector[i] = 0.0f;
    }
}

/**
 * @brief Determines the maximum entry of a vector
 *
 * @param vec Input vector
 * @return float maximum value in vec
 */
float vec_max(const vector_t *const vec)
{
    unsigned short vec_row_size = vec->rows; /*Number of Rows*/

    float value_max = vec->vector[0]; /*Assign first element as start*/

    for (unsigned short i = 1U; i < vec_row_size; i++)
    {
        if (vec->vector[i] > value_max)
        {
            value_max = vec->vector[i];
        }
    }
    return (value_max);
}

/**
 * @brief Vector multiplication with scalar
 *
 * @param vec Vector
 * @param scalar Scalar float
 */
void vec_scalar(vector_t *const vec, float const scalar)
{
    unsigned short vec_row_size = vec->rows; /*Number of Rows*/

    for (unsigned short i = 0U; i < vec_row_size; i++)
    {
        vec->vector[i] = scalar * vec->vector[i];
    }
}

/**
 * @brief Pivoting of a vector with pivot vector
 *
 * @param vec vector to be pivoted (nx1)
 * @param pivot pivoting vector (nx1)
 * @param res_vec result vector (nx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int vec_pivot(const vector_t *const vec, const vector_u_t *const pivot, vector_t *res_vec)
{
    unsigned int status = 0;
    if ((vec->rows == pivot->rows) && (vec->rows == res_vec->rows))
    {
        for (unsigned short i = 0; i < vec->rows; i++)
        {
            res_vec->vector[i] = vec->vector[pivot->vector[i]];
        }
    }
    else
    {
        status = 1U;
    }

    return status;
}

/**
 * @brief Pivoting of a vector with pivot vector
 *
 * @param vec vector to be pivoted (nx1)
 * @param pivot pivoting vector (nx1)
 * @param res_vec result vector (nx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int vec_pivot_c(const vector_c_t *const vec, const vector_u_t *const pivot, vector_t *res_vec)
{
    unsigned int status = 0;
    if ((vec->rows == pivot->rows) && (vec->rows == res_vec->rows))
    {
        for (unsigned short i = 0; i < vec->rows; i++)
        {
            res_vec->vector[i] = vec->vector[pivot->vector[i]];
        }
    }
    else
    {
        status = 1U;
    }

    return status;
}

/**
 * @brief Reverste pivoting of a vector with pivot vector
 *
 * @param vec vector to be pivoted (nx1)
 * @param pivot pivoting vector (nx1)
 * @param res_vec result vector (nx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int vec_pivot_rev(const vector_t *const vec, const vector_u_t *const pivot, vector_t *res_vec)
{
    unsigned int status = 0;
    if ((vec->rows == pivot->rows) && (vec->rows == res_vec->rows))
    {
        for (unsigned short i = 0; i < vec->rows; i++)
        {
            res_vec->vector[pivot->vector[i]] = vec->vector[i];
        }
    }
    else
    {
        status = 1U;
    }

    return status;
}

/**
 * @brief Inverts all elements in a vector
 *
 * @param vec Input vector
 * @return unsigned int 0 = successful exit; 2 = divison by zero
 */
unsigned int vec_inv(vector_t *const vec)
{
    unsigned int status = 0;
    for (unsigned short i = 0; i < vec->rows; i++)
    {
        if (ABS(vec->vector[i]) > EPS_F)
        {
            vec->vector[i] = 1.0f / vec->vector[i];
        }
        else
        {
            status = 2U;
        }
    }
    return (status);
}

/**
 * @brief Calculates the inner product x^T*y
 *
 * @param vec1 vector x (nx1)
 * @param vec2 vector y (nx1)
 * @param res float inner product result
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int vec_T_vec(const vector_t *const vec1, const vector_t *const vec2, float *res)
{
    unsigned int status = 0;
    *res = 0.0f;
    if (vec1->rows != vec2->rows)
    {
        status = 1;
    }
    for (unsigned short i = 0U; i < vec1->rows; i++)
    {
        *res += vec1->vector[i] * vec2->vector[i];
    }
    return (status);
}

/**
 * @brief Calculates the sum of two vectors x + y
 *
 * @param vec1 vector x (nx1)
 * @param vec2 vector y (nx1)
 * @param res_vec vector sum
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int vec_sum(const vector_t *const vec1, const vector_t *const vec2, vector_t *res_vec)
{
    unsigned int status = 0;
    if ((vec1->rows != vec2->rows) && (res_vec->rows != vec2->rows))
    {
        status = 1;
    }
    for (unsigned short i = 0U; i < vec1->rows; i++)
    {
        res_vec->vector[i] = vec1->vector[i] + vec2->vector[i];
    }
    return (status);
}

/**
 * @brief Calculates the sum of two vectors x - y
 *
 * @param vec1 vector x (nx1)
 * @param vec2 vector y (nx1)
 * @param res_vec vector difference
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int vec_diff(const vector_t *const vec1, const vector_t *const vec2, vector_t *res_vec)
{
    unsigned int status = 0;
    if ((vec1->rows != vec2->rows) && (res_vec->rows != vec2->rows))
    {
        status = 1;
    }
    for (unsigned short i = 0U; i < vec1->rows; i++)
    {
        res_vec->vector[i] = vec1->vector[i] - vec2->vector[i];
    }
    return (status);
}

void vec_raw_copy(const float *vec, vector_t *const res_vec)
{
    unsigned short no_of_rows = res_vec->rows; /* Maximum Number of Rows */

    // TODO Check for right dimensions
    for (unsigned short n = 0U; n < no_of_rows; n++)
    {
        res_vec->vector[n] = vec[n];
    }
}

/**
 * @brief Copies over the vector in the result vector
 *
 * @param vec
 * @param res_vec
 */
void vec_copy(const vector_t *const vec, vector_t *const res_vec)
{
    // TODO Check for right dimensions
    for (unsigned short n = 0U; n < vec->rows; n++)
    {
        res_vec->vector[n] = vec->vector[n];
    }
}

void mat_raw_copy(const float *matr, matrix_t *const res_matr)
{
    unsigned short no_of_cols = res_matr->columns; /* Maximum Number of Columns */
    unsigned short no_of_rows = res_matr->rows;    /* Maximum Number of Rows */

    // TODO Check for right dimensions
    for (unsigned short n = 0U; n < no_of_rows; n++)
    {
        for (unsigned short m = 0U; m < no_of_cols; m++)
        {
            M_P_ENTRY(res_matr, n, m) = matr[(m * no_of_rows) + n];
        }
    }
}

/**
 * @brief Extracts the diagonal of a matrix
 *
 * @param matr Input matrix (n x m)
 * @param res_vec Diagonal Vector (smaller of n or m)
 */
void mat_diag(const matrix_t *const matr, vector_t *res_vec)
{
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/

    unsigned short i_max = 0U;
    if (mat_col_size > mat_row_size)
    {
        i_max = mat_row_size;
    }
    else
    {
        i_max = mat_col_size;
    }

    for (unsigned short i = 0U; i < i_max; i++)
    {
        res_vec->vector[i] = M_P_ENTRY(matr, i, i);
    }
}

/**
 * @brief Finds the maximum value on the diagonal of a matrix
 *
 * @param matr Matrix
 * @return float maximum diagonal value
 */
float mat_max_diag(const matrix_t *const matr)
{
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/

    unsigned short i_max = 0U;
    if (mat_col_size > mat_row_size)
    {
        i_max = mat_row_size;
    }
    else
    {
        i_max = mat_col_size;
    }

    float value_max = M_P_ENTRY(matr, 0U, 0U); /*Assign first element as start*/

    for (unsigned short i = 1U; i < i_max; i++)
    {
        if (M_P_ENTRY(matr, i, i) > value_max)
        {
            value_max = M_P_ENTRY(matr, i, i);
        }
    }
    return (value_max);
}

/**
 * @brief Get idx of maximum diagonal value of a matrix matr starting from start
 *
 * @param matr
 * @param start
 * @return unsigned short
 */
unsigned short mat_max_diag_pivot(const matrix_t *const matr, const vector_u_t *const p, unsigned short start)
{
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/

    unsigned short i_max = 0U;
    if (mat_col_size > mat_row_size)
    {
        i_max = mat_row_size;
    }
    else
    {
        i_max = mat_col_size;
    }

    float value_max = M_P_ENTRY(matr, p->vector[start], p->vector[start]); /*Assign first element as start*/
    unsigned short idx_max = start;

    for (unsigned short i = start; i < i_max; i++)
    {
        if (M_P_ENTRY(matr, p->vector[i], p->vector[i]) > value_max)
        {
            value_max = M_P_ENTRY(matr, p->vector[i], p->vector[i]);
            idx_max = i;
        }
    }
    return (idx_max);
}

/**
 * @brief mat_vec implements a matrix vector multiplication $Ax = y$
 *
 * @param matr Matrix A (n x m)
 * @param vec Vector x (m x 1)
 * @param res_vec Vector y (n x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/
    unsigned short vec_row_size = vec->rows;     /* Number of Rows of the Vector */

    if (mat_col_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else // TODO Can be removed
    {
        vec_zero(res_vec);
        for (unsigned short i = 0U; i < mat_row_size; i++)
        {
            for (unsigned short j = 0U; j < mat_col_size; j++)
            {
                res_vec->vector[i] = res_vec->vector[i] + /*Calculate matrix vector product*/
                                     (M_P_ENTRY(matr, i, j) * vec->vector[j]);
            }
        }
        /* return status 0 (OK)         */
    }
    return (status);
}

/**
 * @brief mat_vec implements a matrix transpose vector multiplication $A^Tx = y$
 *
 * @param matr Matrix A (n x m)
 * @param vec Vector x (n x 1)
 * @param res_vec Vector y (m x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_t_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/
    unsigned short vec_row_size = vec->rows;     /* Number of Rows of the Vector */

    if (mat_row_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short j = 0U; j < mat_col_size; j++)
        {
            for (unsigned short i = 0U; i < mat_row_size; i++)
            {
                res_vec->vector[j] = res_vec->vector[j] + /*Calculate matrix vector product*/
                                     (M_P_ENTRY(matr, i, j) * vec->vector[i]);
            }
        }
        /* return status 0 (OK)         */
    }
    return (status);
}

/**
 * @brief mat_vec implements a matrix transpose vector multiplication $A^Tx = y$
 *
 * @param matr const Matrix A (n x m)
 * @param vec Vector x (n x 1)
 * @param res_vec Vector y (m x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_t_vec_c(const matrix_c_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/
    unsigned short vec_row_size = vec->rows;     /* Number of Rows of the Vector */

    if (mat_row_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short j = 0U; j < mat_col_size; j++)
        {
            for (unsigned short i = 0U; i < mat_row_size; i++)
            {
                res_vec->vector[j] = res_vec->vector[j] + /*Calculate matrix vector product*/
                                     (M_P_ENTRY(matr, i, j) * vec->vector[i]);
            }
        }
        /* return status 0 (OK)         */
    }
    return (status);
}

/**
 * @brief Multiplies a diagonal matrix with a vector $D x = y$
 *
 * @param diag Vector (nx1) with the diagonal elements
 * @param vec Vector x (n x 1)
 * @param res_vec Vector y (n x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions;
 */
unsigned int diag_vec(const vector_t *const diag, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short diag_row_size = diag->rows; /* Number of Rows of the Diagonal */
    unsigned short vec_row_size = vec->rows;   /* Number of Rows of the Vector */

    if (diag_row_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        for (unsigned short i = 0U; i < vec_row_size; i++)
        {
            res_vec->vector[i] = vec->vector[i] * diag->vector[i];
        }
        /* return status 0 (OK) */
    }
    return (status);
}

/**
 * @brief Multiplies a diagonal matrix with a vector $D x = y$
 *
 * @param diag const Vector (nx1) with the diagonal elements
 * @param vec const Vector x (n x 1)
 * @param res_vec Vector y (n x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions;
 */
unsigned int diag_vec_c(const vector_c_t *const diag, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short diag_row_size = diag->rows; /* Number of Rows of the Diagonal */
    unsigned short vec_row_size = vec->rows;   /* Number of Rows of the Vector */

    if (diag_row_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        for (unsigned short i = 0U; i < vec_row_size; i++)
        {
            res_vec->vector[i] = vec->vector[i] * diag->vector[i];
        }
        /* return status 0 (OK) */
    }
    return (status);
}

/**
 * @brief Multiplies the inverse of a diagonal matrix with a vector $D^{-1} x = y$
 *
 * @param diag Vector (nx1) with the diagonal elements
 * @param vec Vector x (n x 1)
 * @param res_vec Vector y (n x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions; 2 = divison by zero
 */
unsigned int diag_inv_vec(const vector_t *const diag, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short diag_row_size = diag->rows; /* Number of Rows of the Diagonal */
    unsigned short vec_row_size = vec->rows;   /* Number of Rows of the Vector */

    if (diag_row_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short i = 0U; i < vec_row_size; i++)
        {
            if (ABS(diag->vector[i]) > EPS_F)
            {
                res_vec->vector[i] = vec->vector[i] / diag->vector[i];
            }
            else
            {
                status = 2U; /* Divison by 0 */
                break;
            }
        }
    }
    return (status); /* return status 0 (OK) */
}

/**
 * @brief lower_vec implements a lower matrix vector multiplication $Lx = y$. Also working for non squared lower matrices
 *
 * @param matr Matrix L (n x m) which is a lower matrix
 * @param vec Vector x (m x 1)
 * @param res_vec Vector y (n x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int lower_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/
    unsigned short vec_row_size = vec->rows;     /* Number of Rows of the Vector */

    if (mat_col_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short i = 0U; i < mat_row_size; i++)
        {
            unsigned short max_j = 0U;
            /*Min of i and mat_col_size*/
            if (i >= mat_col_size)
            {
                max_j = mat_col_size;
            }
            else
            {
                max_j = i + 1U;
            }

            /*Loop through lower*/
            for (unsigned short j = 0U; j < max_j; j++)
            {
                res_vec->vector[i] = res_vec->vector[i] + /*Calculate matrix vector product*/
                                     (M_P_ENTRY(matr, i, j) * vec->vector[j]);
            }
        }
        /* return status 0 (OK) */
    }
    return (status);
}

/**
 * @brief lower_t_vec implements a lower matrix transpose vector multiplication $L^Tx = y$. Also working for non squared lower matrices
 *
 * @param matr Matrix L (n x m) which is a lower matrix
 * @param vec Vector x (n x 1)
 * @param res_vec Vector y (m x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int lower_t_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/
    unsigned short vec_row_size = vec->rows;     /* Number of Rows of the Vector */

    if (mat_row_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        unsigned short max_j = 0U;

        /*Min of mat_row_size and mat_col_size*/
        if (mat_col_size > mat_row_size)
        {
            max_j = mat_row_size;
        }
        else
        {
            max_j = mat_col_size;
        }

        for (unsigned short j = 0U; j < max_j; j++)
        {
            for (unsigned short i = j; i < mat_row_size; i++)
            {
                res_vec->vector[j] = res_vec->vector[j] + /*Calculate matrix vector product*/
                                     (M_P_ENTRY(matr, i, j) * vec->vector[i]);
            }
        }
        /* return status 0 (OK)         */
    }
    return (status);
}

/**
 * @brief lower_t_vec implements a lower matrix transpose vector multiplication $L^Tx = y$. Also working for non squared lower matrices
 *
 * @param matr const Matrix L (n x m) which is a lower matrix
 * @param vec const Vector x (n x 1)
 * @param res_vec Vector y (m x 1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int lower_t_vec_c(const matrix_t *const matr, const vector_c_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/
    unsigned short vec_row_size = vec->rows;     /* Number of Rows of the Vector */

    if (mat_row_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        unsigned short max_j = 0U;

        /*Min of mat_row_size and mat_col_size*/
        if (mat_col_size > mat_row_size)
        {
            max_j = mat_row_size;
        }
        else
        {
            max_j = mat_col_size;
        }

        for (unsigned short j = 0U; j < max_j; j++)
        {
            for (unsigned short i = j; i < mat_row_size; i++)
            {
                res_vec->vector[j] = res_vec->vector[j] + /*Calculate matrix vector product*/
                                     (M_P_ENTRY(matr, i, j) * vec->vector[i]);
            }
        }
        /* return status 0 (OK)         */
    }
    return (status);
}

/**
 * @brief Calculates the inverse of a lower matrix explicitly
 *
 * @param matr Lower Matrix to be inverted
 * @param res_matr Inverse of matr
 * @return unsigned int unsigned int 0 = successful exit; 1 = wrong dimensions; 2 matrix not invertible
 */
unsigned int lower_inv(const matrix_t *const matr, matrix_t *res_matr)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/

    if (mat_col_size != mat_row_size)
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        eye(res_matr, mat_col_size);

        for (unsigned short i = 0U; i < mat_row_size; i++)
        {
            if ((M_P_ENTRY(matr, i, i) < EPS_F) && (M_P_ENTRY(matr, i, i) > -EPS_F))
            {
                status = 2U;
                break;
            }

            M_P_ENTRY(res_matr, i, i) = 1.0f / M_P_ENTRY(matr, i, i);

            /*Normalizing the row by the diag element*/
            for (unsigned short j = 0U; j < i; j++)
            {
                M_P_ENTRY(matr, i, j) = M_P_ENTRY(matr, i, j) / M_P_ENTRY(matr, i, i);
            }

            /*Bringing the elements left to the diagonal to zero*/
            if (i > 0U) /*Only relevant from second row*/
            {
                for (unsigned short idx = 0U; idx < i; idx++)
                {
                    for (unsigned short j = 0U; j < (idx + 1U); j++)
                    {
                        M_P_ENTRY(res_matr, i, j) = M_P_ENTRY(res_matr, i, j) -
                                                    (M_P_ENTRY(res_matr, idx, j) * M_P_ENTRY(matr, i, idx));
                    }
                }
            }
        }
    }
    return (status);
}

/**
 * @brief Calculates the inverse of a unit lower matrix explicitly
 *
 * @param matr Unit Lower Matrix to be inverted
 * @param res_matr Inverse of matr
 * @return unsigned int unsigned int 0 = successful exit; 1 = wrong dimensions; 2 matrix not invertible
 */
unsigned int lower_inv_LDLT(const matrix_t *const matr, matrix_t *res_matr)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/

    if (mat_col_size != mat_row_size)
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        eye(res_matr, mat_col_size);

        for (unsigned short i = 0U; i < mat_row_size; i++)
        {
            if ((M_P_ENTRY(matr, i, i) < EPS_F) && (M_P_ENTRY(matr, i, i) > -EPS_F))
            {
                status = 2U;
                break;
            }

            /*Bringing the elements left to the diagonal to zero*/
            if (i > 0U) /*Only relevant from second row*/
            {
                for (unsigned short idx = 0U; idx < i; idx++)
                {
                    for (unsigned short j = 0U; j < (idx + 1U); j++)
                    {
                        M_P_ENTRY(res_matr, i, j) = M_P_ENTRY(res_matr, i, j) -
                                                    (M_P_ENTRY(res_matr, idx, j) * M_P_ENTRY(matr, i, idx));
                    }
                }
            }
        }
    }
    return (status);
}

/**
 * @brief Calculates the determinant of a 2x2 submatrix
 *
 * @param matr Complete matrix
 * @param row_1 Row index 1 of the submatrix
 * @param row_2 Row index 2 of the submatrix
 * @param col_1 Column index 1 of the submatrix
 * @param col_2 Column index 2 of the submatrix
 * @return float float Determinant of the 2x2 submatrix
 */
float det_two(const matrix_t *const matr, const unsigned short row_1, const unsigned short row_2, const unsigned short col_1, const unsigned short col_2)
{
    float det = 0.0f;

    det = (M_P_ENTRY(matr, row_1, col_1) * M_P_ENTRY(matr, row_2, col_2)) - (M_P_ENTRY(matr, row_1, col_2) * M_P_ENTRY(matr, row_2, col_1));

    return (det);
}

/**
 * @brief Calculates the determinant of a 3x3 submatrix
 *
 * @param matr Complete matrix
 * @param row_1 Row index 1 of the submatrix
 * @param row_2 Row index 2 of the submatrix
 * @param row_3 Row index 3 of the submatrix
 * @param col_1 Column index 1 of the submatrix
 * @param col_2 Column index 2 of the submatrix
 * @param col_3 Column index 3 of the submatrix
 * @return float float Determinant of the 3x3 submatrix
 */
float det_three(const matrix_t *const matr, const unsigned short row_1, const unsigned short row_2, const unsigned short row_3, const unsigned short col_1, const unsigned short col_2, const unsigned short col_3)
{
    float det = 0.0f;

    det = (M_P_ENTRY(matr, row_1, col_1) *
           ((M_P_ENTRY(matr, row_2, col_2) * M_P_ENTRY(matr, row_3, col_3)) - (M_P_ENTRY(matr, row_2, col_3) * M_P_ENTRY(matr, row_3, col_2)))) +
          (M_P_ENTRY(matr, row_1, col_2) *
           ((M_P_ENTRY(matr, row_2, col_3) * M_P_ENTRY(matr, row_3, col_1)) - (M_P_ENTRY(matr, row_2, col_1) * M_P_ENTRY(matr, row_3, col_3)))) +
          (M_P_ENTRY(matr, row_1, col_3) *
           ((M_P_ENTRY(matr, row_2, col_1) * M_P_ENTRY(matr, row_3, col_2)) - (M_P_ENTRY(matr, row_2, col_2) * M_P_ENTRY(matr, row_3, col_1))));

    return (det);
}

/**
 * @brief Implements B*B^T
 *
 * @param matr (n x m)
 * @param W (m x 1)
 * @param res_matr (n x n)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_gram(const matrix_t *const matr, matrix_t *res_matr)
{
    unsigned int status = 0;
    unsigned short matr_col_size; /* Number of Columns of Matrix */
    unsigned short matr_row_size; /* Number of Rows of Matrix */

    matr_col_size = matr->columns; /* get size of matrix */
    matr_row_size = matr->rows;

    mat_zero(res_matr);
    for (unsigned short h = 0U; h < matr_row_size; h++) /* left row counter */
    {
        for (unsigned short j = 0U; j < matr_row_size; j++) /*  right column counter */
        {
            if (j >= h)
            {
                for (unsigned short i = 0U; i < matr_col_size; i++) /* left col counter */
                {
                    M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, h, j) + /* calculate result */
                                                (M_P_ENTRY(matr, h, i) * M_P_ENTRY(matr, j, i));
                }
            }
            else
            {
                M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, j, h);
            }
        }
    }
    /* return status 0 (OK) */
    return (status);
}

/**
 * @brief Implements B*diag(W)*B^T
 *
 * @param matr (n x m)
 * @param W (m x 1)
 * @param res_matr (n x n)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_gram_weighted(const matrix_t *const matr, const vector_t *const W, matrix_t *res_matr)
{
    unsigned int status = 0;
    unsigned short matr_col_size; /* Number of Columns of Matrix */
    unsigned short matr_row_size; /* Number of Rows of Matrix */
    unsigned short vec_row_size;  /* Number of Rows of Vector */

    matr_col_size = matr->columns; /* get size of matrix */
    matr_row_size = matr->rows;

    vec_row_size = W->rows; /* get size of vector */

    if (matr_col_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        mat_zero(res_matr);
        for (unsigned short h = 0U; h < matr_row_size; h++) /* left row counter */
        {
            for (unsigned short j = 0U; j < matr_row_size; j++) /*  right column counter */
            {
                if (j >= h)
                {
                    for (unsigned short i = 0U; i < matr_col_size; i++) /* left col counter */
                    {
                        M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, h, j) + /* calculate result */
                                                    (M_P_ENTRY(matr, h, i) * W->vector[i] * M_P_ENTRY(matr, j, i));
                    }
                }
                else
                {
                    M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, j, h);
                }
            }
        }
        /* return status 0 (OK) */
    }
    return (status);
}

/**
 * @brief Implements B*diag(W)*B^T with constant B and W
 *
 * @param matr (n x m)
 * @param W (m x 1)
 * @param res_matr (n x n)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_gram_weighted_c(const matrix_c_t *const matr, const vector_c_t *const W, matrix_t *res_matr)
{
    unsigned int status = 0;
    unsigned short matr_col_size; /* Number of Columns of Matrix */
    unsigned short matr_row_size; /* Number of Rows of Matrix */
    unsigned short vec_row_size;  /* Number of Rows of Vector */

    matr_col_size = matr->columns; /* get size of matrix */
    matr_row_size = matr->rows;

    vec_row_size = W->rows; /* get size of vector */

    if (matr_col_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        mat_zero(res_matr);
        for (unsigned short h = 0U; h < matr_row_size; h++) /* left row counter */
        {
            for (unsigned short j = 0U; j < matr_row_size; j++) /*  right column counter */
            {
                if (j >= h)
                {
                    for (unsigned short i = 0U; i < matr_col_size; i++) /* left col counter */
                    {
                        M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, h, j) + /* calculate result */
                                                    (M_P_ENTRY(matr, h, i) * W->vector[i] * M_P_ENTRY(matr, j, i));
                    }
                }
                else
                {
                    M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, j, h);
                }
            }
        }
        /* return status 0 (OK) */
    }
    return (status);
}

/**
 * @brief Implements B^T*diag(W)*B
 *
 * @param matr (n x m)
 * @param W (n x 1)
 * @param res_matr (n x n)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_gram_T_weighted(const matrix_t *const matr, const vector_t *const W, matrix_t *res_matr)
{
    unsigned int status = 0;
    unsigned short matr_col_size; /* Number of Columns of Matrix */
    unsigned short matr_row_size; /* Number of Rows of Matrix */
    unsigned short vec_row_size;  /* Number of Rows of Vector */

    matr_col_size = matr->columns; /* get size of matrix */
    matr_row_size = matr->rows;

    vec_row_size = W->rows; /* get size of vector */

    if (matr_row_size != vec_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        mat_zero(res_matr);
        for (unsigned short h = 0U; h < matr_col_size; h++) /* left row counter */
        {
            for (unsigned short j = 0U; j < matr_col_size; j++) /*  right column counter */
            {
                if (j >= h)
                {
                    for (unsigned short i = 0U; i < matr_row_size; i++) /* left col counter */
                    {
                        M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, h, j) + /* calculate result */
                                                    (M_P_ENTRY(matr, i, h) * W->vector[i] * M_P_ENTRY(matr, i, j));
                    }
                }
                else
                {
                    M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, j, h);
                }
            }
        }
        /* return status 0 (OK) */
    }
    return (status);
}

/**
 * @brief Implements a full matrix times a diagonal matrix matr*diag(W)
 *
 * @param matr (nxm)
 * @param W (mx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int mat_mult_diag_c(matrix_t *const matr, const vector_c_t *const W)
{
    unsigned int status = 0;
    unsigned short mat_col_size = matr->columns; /* Number of Columns of the Matrix*/
    unsigned short mat_row_size = matr->rows;    /* Number of Rows of the Matrix*/
    unsigned short diag_row_size = W->rows;      /* Number of Rows of the Vector */

    if (mat_col_size != diag_row_size) /*  check if sizes do match... */
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else // TODO Can be removed
    {
        for (unsigned short i = 0U; i < mat_row_size; i++)
        {
            for (unsigned short j = 0U; j < mat_col_size; j++)
            {
                M_P_ENTRY(matr, i, j) = M_P_ENTRY(matr, i, j) * W->vector[j]; /* Scale each column by the corresponding element in W*/
            }
        }
        /* return status 0 (OK)  */
    }
    return (status);
}

/**
 * @brief Implements L^T*L also for lower L with more rows the columns
 *
 * @param matr Lower matrix L(nxm) can have more rows the columns
 * @param res_matr mxm
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int lower_T_lower(const matrix_t *const matr, matrix_t *res_matr)
{
    unsigned int status = 0;
    unsigned short matr_col_size; /* Number of Columns of Matrix */
    unsigned short matr_row_size; /* Number of Rows of Matrix */

    matr_col_size = matr->columns; /* get size of matrix */
    matr_row_size = matr->rows;
    if (matr_col_size > matr_row_size)
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        mat_zero(res_matr);
        for (unsigned short h = 0U; h < matr_col_size; h++) /* left row counter */
        {
            for (unsigned short j = 0U; j < matr_col_size; j++) /*  right column counter */
            {
                if (j >= h)
                {
                    for (unsigned short i = j; i < matr_row_size; i++) /* left col counter */
                    {
                        M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, h, j) + /* calculate result */
                                                    (M_P_ENTRY(matr, i, h) * M_P_ENTRY(matr, i, j));
                    }
                }
                else
                {
                    M_P_ENTRY(res_matr, h, j) = M_P_ENTRY(res_matr, j, h);
                }
            }
        }
    }
    return (status); /* return status 0 (OK) */
}

/**
 * @brief Solves the equation L*x=y by forward substitution
 *
 * @param matr Lower matrix L (nxn)
 * @param vec vector y (nx1)
 * @param res_vec vector x (nx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions; 2 = rank deficient
 */
unsigned int solve_lower_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short matr_row_size = matr->rows;
    if (matr_row_size != vec->rows)
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short i = 0U; i < matr_row_size; i++)
        {
            for (unsigned short j = 0U; j < i; j++)
            {
                res_vec->vector[i] += M_P_ENTRY(matr, i, j) * res_vec->vector[j];
            }
            if (ABS(M_P_ENTRY(matr, i, i)) > EPS_F)
            {
                res_vec->vector[i] = (vec->vector[i] - res_vec->vector[i]) / M_P_ENTRY(matr, i, i);
            }
            else
            {
                status = 2U;
                break;
            }
        }
    }
    return (status);
}

/**
 * @brief Solves the equation L^T*x=y by forward substitution
 *
 * @param matr Lower matrix L (nxn)
 * @param vec vector y (nx1)
 * @param res_vec vector x (nx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions; 2 = rank deficient
 */
unsigned int solve_lower_T_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short matr_row_size = matr->rows;
    if (matr_row_size != vec->rows)
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short i = (matr_row_size - 1U); i < matr_row_size; i--)
        {
            for (unsigned short j = (matr_row_size - 1U); j > i; j--)
            {
                res_vec->vector[i] += M_P_ENTRY(matr, j, i) * res_vec->vector[j];
            }
            if (ABS(M_P_ENTRY(matr, i, i)) > EPS_F)
            {
                res_vec->vector[i] = (vec->vector[i] - res_vec->vector[i]) / M_P_ENTRY(matr, i, i);
            }
            else
            {
                status = 2U;
                break;
            }
        }
    }
    return (status);
}

/**
 * @brief Solves the equation L*x=y of a LDL^T decomposition (ones on diag) by forward substitution
 *
 * @param matr Lower matrix L (nxn)
 * @param vec vector y (nx1)
 * @param res_vec vector x (nx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions; 2 = rank deficient
 */
unsigned int solve_LDL_lower_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short matr_row_size = matr->rows;
    if (matr_row_size != vec->rows)
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short i = 0U; i < matr_row_size; i++)
        {
            for (unsigned short j = 0U; j < i; j++)
            {
                res_vec->vector[i] += M_P_ENTRY(matr, i, j) * res_vec->vector[j];
            }
            if (ABS(M_P_ENTRY(matr, i, i)) < EPS_F)
            {
                status = 2U;
                break;
            }
            else
            {
                res_vec->vector[i] = vec->vector[i] - res_vec->vector[i];
            }
        }
    }
    return (status);
}

/**
 * @brief Solves the equation L*x=y of a LDL^T decomposition (ones on diag) by forward substitution
 *
 * @param matr Lower matrix L (nxn)
 * @param vec constant vector y (nx1)
 * @param res_vec vector x (nx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions; 2 = rank deficient
 */
unsigned int solve_LDL_lower_vec_c(const matrix_t *const matr, const vector_c_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short matr_row_size = matr->rows;
    if (matr_row_size != vec->rows)
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short i = 0U; i < matr_row_size; i++)
        {
            for (unsigned short j = 0U; j < i; j++)
            {
                res_vec->vector[i] += M_P_ENTRY(matr, i, j) * res_vec->vector[j];
            }
            if (ABS(M_P_ENTRY(matr, i, i)) < EPS_F)
            {
                status = 2U;
                break;
            }
            else
            {
                res_vec->vector[i] = vec->vector[i] - res_vec->vector[i];
            }
        }
    }
    return (status);
}

/**
 * @brief Solves the equation L^T*x=y of a LDL^T decomposition (ones on diag) by forward substitution
 *
 * @param matr Lower matrix L (nxn)
 * @param vec vector y (nx1)
 * @param res_vec vector x (nx1)
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions; 2 = rank deficient
 */
unsigned int solve_LDL_lower_T_vec(const matrix_t *const matr, const vector_t *const vec, vector_t *res_vec)
{
    unsigned int status = 0;
    unsigned short matr_row_size = matr->rows;
    if (matr_row_size != vec->rows)
    {
        status = 1; /* return status 1 (ERROR) */
    }
    else
    {
        vec_zero(res_vec);
        for (unsigned short i = (matr_row_size - 1U); i < matr_row_size; i--)
        {
            for (unsigned short j = (matr_row_size - 1U); j > i; j--)
            {
                res_vec->vector[i] += M_P_ENTRY(matr, j, i) * res_vec->vector[j];
            }
            if (ABS(M_P_ENTRY(matr, i, i)) < EPS_F)
            {
                status = 2U;
                break;
            }
            else
            {
                res_vec->vector[i] = vec->vector[i] - res_vec->vector[i];
            }
        }
    }
    return (status);
}