/**
 * Copyright (C) 2023-2024 Simon Hafner - All Rights Reserved  
 * Institute of Flight System Dynamics -- Technische Universitaet Muenchen
 * Contact: simon.hafner@tum.de  
 */

#include "sym_mat_inv.h"

/**
 * @brief Calculates the inverse of a general 2x2 matrix based on the adjoint and determinat
 *
 * @param matr Input matrix (2x2)
 * @param res_matr Inverse matrix (2x2)
 * @return unsigned int 0U = successful exit; 2U = not invertible
 */
unsigned int sym_inv_two(const matrix_t *const matr, matrix_t *res_matr)
{
    unsigned int status = 0U;
    float determinant;

    // TODO Add check of input matrix

    res_matr->matrix[0U] = +matr->matrix[matr->rows + 1U];
    res_matr->matrix[res_matr->rows] = -matr->matrix[1U];

    res_matr->matrix[1U] = -matr->matrix[matr->rows];
    res_matr->matrix[res_matr->rows + 1U] = +matr->matrix[0U];

    determinant = (matr->matrix[0U] * res_matr->matrix[0U]) +
                  +(matr->matrix[1U] * res_matr->matrix[res_matr->rows]);

    if (ABS(determinant) < EPS_F)
    {
        status = 2U;
    }
    else
    {
        mat_mult_scalar(res_matr, 1.0f / determinant);
    }
    return (status);
}

/**
 * @brief Calculates the inverse of a symmetric 3x3 matrix based on the adjoint and determinat
 *
 * @param matr Input matrix (3x3)
 * @param res_matr Inverse matrix (3x3)
 * @return unsigned int 0U = successful exit; 2U = not invertible
 */
unsigned int sym_inv_three(const matrix_t *const matr, matrix_t *res_matr)
{
    unsigned int status = 0U;
    float determinant;

    // TODO Add check of input matrix
    res_matr->matrix[0U] = +det_two(matr, 1U, 2U, 1U, 2U);
    res_matr->matrix[res_matr->rows] = -det_two(matr, 1U, 2U, 0U, 2U);
    res_matr->matrix[2U * res_matr->rows] = +det_two(matr, 1U, 2U, 0U, 1U);

    res_matr->matrix[1U] = res_matr->matrix[res_matr->rows];
    res_matr->matrix[res_matr->rows + 1U] = +det_two(matr, 0U, 2U, 0U, 2U);
    res_matr->matrix[(2U * res_matr->rows) + 1U] = -det_two(matr, 0U, 2U, 0U, 1U);

    res_matr->matrix[2U] = res_matr->matrix[2U * res_matr->rows];
    res_matr->matrix[res_matr->rows + 2U] = res_matr->matrix[(2U * res_matr->rows) + 1U];
    res_matr->matrix[(2U * res_matr->rows) + 2U] = +det_two(matr, 0U, 1U, 0U, 1U);

    // Determinant
    determinant = (matr->matrix[0U] * res_matr->matrix[0U]) + (matr->matrix[matr->rows] * res_matr->matrix[res_matr->rows]) + (matr->matrix[2U * matr->rows] * res_matr->matrix[2U * res_matr->rows]);

    if (ABS(determinant) < EPS_F)
    {
        status = 2U;
    }
    else
    {
        mat_mult_scalar(res_matr, 1.0f / determinant);
    }
    return (status);
}

/**
 * @brief Calculates the inverse of a symmetric 4x4 matrix based on the adjoint and determinat
 *
 * @param matr Input matrix (4x4)
 * @param res_matr Inverse matrix (4x4)
 * @return unsigned int 0U = successful exit; 1U = not invertible
 */
unsigned int sym_inv_four(const matrix_t *const matr, matrix_t *res_matr)
{
    unsigned int status = 0U;
    float determinant;

    // TODO Add check of input matrix

    res_matr->matrix[0U] = +det_three(matr, 1U, 2U, 3U, 1U, 2U, 3U);
    res_matr->matrix[res_matr->rows] = -det_three(matr, 1U, 2U, 3U, 0U, 2U, 3U);
    res_matr->matrix[2U * res_matr->rows] = +det_three(matr, 1U, 2U, 3U, 0U, 1U, 3U);
    res_matr->matrix[3U * res_matr->rows] = -det_three(matr, 1U, 2U, 3U, 0U, 1U, 2U);

    res_matr->matrix[1U] = res_matr->matrix[res_matr->rows];
    res_matr->matrix[res_matr->rows + 1U] = +det_three(matr, 0U, 2U, 3U, 0U, 2U, 3U);
    res_matr->matrix[(2U * res_matr->rows) + 1U] = -det_three(matr, 0U, 2U, 3U, 0U, 1U, 3U);
    res_matr->matrix[(3U * res_matr->rows) + 1U] = +det_three(matr, 0U, 2U, 3U, 0U, 1U, 2U);

    res_matr->matrix[2U] = res_matr->matrix[2U * res_matr->rows];
    res_matr->matrix[res_matr->rows + 2U] = res_matr->matrix[(2U * res_matr->rows) + 1U];
    res_matr->matrix[(2U * res_matr->rows) + 2U] = +det_three(matr, 0U, 1U, 3U, 0U, 1U, 3U);
    res_matr->matrix[(3U * res_matr->rows) + 2U] = -det_three(matr, 0U, 1U, 3U, 0U, 1U, 2U);

    res_matr->matrix[3U] = res_matr->matrix[3U * res_matr->rows];
    res_matr->matrix[res_matr->rows + 3U] = res_matr->matrix[(3U * res_matr->rows) + 1U];
    res_matr->matrix[(2U * res_matr->rows) + 3U] = res_matr->matrix[(3U * res_matr->rows) + 2U];
    res_matr->matrix[(3U * res_matr->rows) + 3U] = +det_three(matr, 0U, 1U, 2U, 0U, 1U, 2U);

    /* Determinant */
    determinant = (matr->matrix[0U] * res_matr->matrix[0U]) + (matr->matrix[matr->rows] * res_matr->matrix[res_matr->rows]) +
                  (matr->matrix[2U * matr->rows] * res_matr->matrix[2U * res_matr->rows]) + (matr->matrix[3U * matr->rows] * res_matr->matrix[3U * res_matr->rows]);

    if (ABS(determinant) < EPS_F)
    {
        status = 2U;
    }
    else
    {
        mat_mult_scalar(res_matr, 1.0f / determinant);
    }
    return (status);
}

/**
 * @brief Calculates the inverse of a symmetric 5x5 matrix based on the adjoint and determinat
 *
 * @param matr Input matrix (5x5)
 * @param res_matr Inverse matrix (5x5)
 * @return unsigned int 0U = successful exit; 1U = not invertible
 */
unsigned int sym_inv_five(const matrix_t *const matr, matrix_t *res_matr)
{
    unsigned int status = 0U;
    float determinant;

    // TODO Add check of input matrix

    float det_3_1 = det_three(matr, 2U, 3U, 4U, 2U, 3U, 4U);
    float det_3_2 = det_three(matr, 2U, 3U, 4U, 1U, 3U, 4U);
    float det_3_3 = det_three(matr, 2U, 3U, 4U, 1U, 2U, 4U);
    float det_3_4 = det_three(matr, 2U, 3U, 4U, 1U, 2U, 3U);
    float det_3_6 = det_three(matr, 2U, 3U, 4U, 0U, 3U, 4U);
    float det_3_7 = det_three(matr, 2U, 3U, 4U, 0U, 2U, 4U);
    float det_3_8 = det_three(matr, 2U, 3U, 4U, 0U, 2U, 3U);
    float det_3_11 = det_three(matr, 2U, 3U, 4U, 0U, 1U, 4U);
    float det_3_12 = det_three(matr, 2U, 3U, 4U, 0U, 1U, 3U);
    float det_3_16 = det_three(matr, 2U, 3U, 4U, 0U, 1U, 2U);
    float det_3_37 = det_three(matr, 1U, 3U, 4U, 1U, 3U, 4U);
    float det_3_38 = det_three(matr, 1U, 3U, 4U, 0U, 3U, 4U);
    float det_3_39 = det_three(matr, 1U, 3U, 4U, 0U, 1U, 4U);
    float det_3_40 = det_three(matr, 1U, 3U, 4U, 0U, 1U, 3U);
    float det_3_41 = det_three(matr, 1U, 3U, 4U, 1U, 2U, 4U);
    float det_3_42 = det_three(matr, 1U, 3U, 4U, 0U, 2U, 4U);
    float det_3_44 = det_three(matr, 1U, 3U, 4U, 0U, 1U, 2U);
    float det_3_45 = det_three(matr, 1U, 3U, 4U, 1U, 2U, 3U);
    float det_3_46 = det_three(matr, 1U, 3U, 4U, 0U, 2U, 3U);
    float det_3_49 = det_three(matr, 1U, 2U, 4U, 1U, 2U, 4U);
    float det_3_50 = det_three(matr, 1U, 2U, 4U, 0U, 2U, 4U);
    float det_3_51 = det_three(matr, 1U, 2U, 4U, 0U, 1U, 4U);
    float det_3_52 = det_three(matr, 1U, 2U, 4U, 0U, 1U, 2U);
    float det_3_53 = det_three(matr, 1U, 2U, 4U, 1U, 2U, 3U);
    float det_3_54 = det_three(matr, 1U, 2U, 4U, 0U, 2U, 3U);
    float det_3_55 = det_three(matr, 1U, 2U, 4U, 0U, 1U, 3U);
    float det_3_57 = det_three(matr, 1U, 2U, 3U, 1U, 2U, 3U);
    float det_3_58 = det_three(matr, 1U, 2U, 3U, 0U, 2U, 3U);
    float det_3_59 = det_three(matr, 1U, 2U, 3U, 0U, 1U, 3U);
    float det_3_60 = det_three(matr, 1U, 2U, 3U, 0U, 1U, 2U);

    // Adjoint

    res_matr->matrix[0U] = +(+(matr->matrix[matr->rows + 1U] * det_3_1) - (matr->matrix[(2U * matr->rows) + 1U] * det_3_2) + (matr->matrix[(3U * matr->rows) + 1U] * det_3_3) - (matr->matrix[(4U * matr->rows) + 1U] * det_3_4));
    res_matr->matrix[res_matr->rows] = -(+(matr->matrix[1U] * det_3_1) - (matr->matrix[(2U * matr->rows) + 1U] * det_3_6) + (matr->matrix[(3U * matr->rows) + 1U] * det_3_7) - (matr->matrix[(4U * matr->rows) + 1U] * det_3_8));
    res_matr->matrix[2U * res_matr->rows] = +(+(matr->matrix[1U] * det_3_2) - (matr->matrix[matr->rows + 1U] * det_3_6) + (matr->matrix[(3U * matr->rows) + 1U] * det_3_11) - (matr->matrix[(4U * matr->rows) + 1U] * det_3_12));
    res_matr->matrix[3U * res_matr->rows] = -(+(matr->matrix[1U] * det_3_3) - (matr->matrix[matr->rows + 1U] * det_3_7) + (matr->matrix[(2U * matr->rows) + 1U] * det_3_11) - (matr->matrix[(4U * matr->rows) + 1U] * det_3_16));
    res_matr->matrix[4U * res_matr->rows] = +(+(matr->matrix[1U] * det_3_4) - (matr->matrix[matr->rows + 1U] * det_3_8) + (matr->matrix[(2U * matr->rows) + 1U] * det_3_12) - (matr->matrix[(3U * matr->rows) + 1U] * det_3_16));

    res_matr->matrix[1U] = res_matr->matrix[res_matr->rows];
    res_matr->matrix[res_matr->rows + 1U] = +(+(matr->matrix[0U] * det_3_1) - (matr->matrix[2U * matr->rows] * det_3_6) + (matr->matrix[3U * matr->rows] * det_3_7) - (matr->matrix[4U * matr->rows] * det_3_8));
    res_matr->matrix[(2U * res_matr->rows) + 1U] = -(+(matr->matrix[0U] * det_3_2) - (matr->matrix[matr->rows] * det_3_6) + (matr->matrix[3U * matr->rows] * det_3_11) - (matr->matrix[4U * matr->rows] * det_3_12));
    res_matr->matrix[(3U * res_matr->rows) + 1U] = +(+(matr->matrix[0U] * det_3_3) - (matr->matrix[matr->rows] * det_3_7) + (matr->matrix[2U * matr->rows] * det_3_11) - (matr->matrix[4U * matr->rows] * det_3_16));
    res_matr->matrix[(4U * res_matr->rows) + 1U] = -(+(matr->matrix[0U] * det_3_4) - (matr->matrix[matr->rows] * det_3_8) + (matr->matrix[2U * matr->rows] * det_3_12) - (matr->matrix[3U * matr->rows] * det_3_16));

    res_matr->matrix[2U] = res_matr->matrix[2U * res_matr->rows];
    res_matr->matrix[res_matr->rows + 2U] = res_matr->matrix[(2U * res_matr->rows) + 1U];
    res_matr->matrix[(2U * res_matr->rows) + 2U] = +(+(matr->matrix[0U] * det_3_37) - (matr->matrix[matr->rows] * det_3_38) + (matr->matrix[3U * matr->rows] * det_3_39) - (matr->matrix[4U * matr->rows] * det_3_40));
    res_matr->matrix[(3U * res_matr->rows) + 2U] = -(+(matr->matrix[0U] * det_3_41) - (matr->matrix[matr->rows] * det_3_42) + (matr->matrix[2U * matr->rows] * det_3_39) - (matr->matrix[4U * matr->rows] * det_3_44));
    res_matr->matrix[(4U * res_matr->rows) + 2U] = +(+(matr->matrix[0U] * det_3_45) - (matr->matrix[matr->rows] * det_3_46) + (matr->matrix[2U * matr->rows] * det_3_40) - (matr->matrix[3U * matr->rows] * det_3_44));

    res_matr->matrix[3U] = res_matr->matrix[3U * res_matr->rows];
    res_matr->matrix[res_matr->rows + 3U] = res_matr->matrix[(3U * res_matr->rows) + 1U];
    res_matr->matrix[(2U * res_matr->rows) + 3U] = res_matr->matrix[(3U * res_matr->rows) + 2U];
    res_matr->matrix[(3U * res_matr->rows) + 3U] = +(+(matr->matrix[0U] * det_3_49) - (matr->matrix[matr->rows] * det_3_50) + (matr->matrix[2U * matr->rows] * det_3_51) - (matr->matrix[4U * matr->rows] * det_3_52));
    res_matr->matrix[(4U * res_matr->rows) + 3U] = -(+(matr->matrix[0U] * det_3_53) - (matr->matrix[matr->rows] * det_3_54) + (matr->matrix[2U * matr->rows] * det_3_55) - (matr->matrix[3U * matr->rows] * det_3_52));

    res_matr->matrix[4U] = res_matr->matrix[4U * res_matr->rows];
    res_matr->matrix[res_matr->rows + 4U] = res_matr->matrix[(4U * res_matr->rows) + 1U];
    res_matr->matrix[(2U * res_matr->rows) + 4U] = res_matr->matrix[(4U * res_matr->rows) + 2U];
    res_matr->matrix[(3U * res_matr->rows) + 4U] = res_matr->matrix[(4U * res_matr->rows) + 3U];
    res_matr->matrix[(4U * res_matr->rows) + 4U] = +(+(matr->matrix[0U] * det_3_57) - (matr->matrix[matr->rows] * det_3_58) + (matr->matrix[2U * matr->rows] * det_3_59) - (matr->matrix[3U * matr->rows] * det_3_60));

    // Determinant
    determinant = (matr->matrix[0U] * res_matr->matrix[0U]) + (matr->matrix[matr->rows] * res_matr->matrix[res_matr->rows]) +
                  (matr->matrix[2U * matr->rows] * res_matr->matrix[2U * res_matr->rows]) + (matr->matrix[3U * matr->rows] * res_matr->matrix[3U * res_matr->rows]) +
                  (matr->matrix[4U * matr->rows] * res_matr->matrix[4U * res_matr->rows]);

    if (ABS(determinant) < EPS_F)
    {
        status = 2U;
    }
    else
    {
        mat_mult_scalar(res_matr, 1.0f / determinant);
    }
    return (status);
}
