/*------------------------------------------------------------------------*/
/*                       control_allocation_LDL.c                         */
/*------------------------------------------------------------------------*/
/*                                                                        */
/*  I N S T I T U T E   O F   F L I G H T   S Y S T E M   D Y N A M I C S */
/*                       Web: www.fsd.ed.tum.de                           */
/*                         _______                                        */
/*                            |  |    |\  /|                              */
/*                            |  |    | \/ |                              */
/*                            |  |____|    |                              */
/*                 Technische Universitaet Muenchen TUM                   */
/*                                                                        */
/*(c) 2023 by Institute of Flight System Dynamics                         */
/*                        All Rights Reserved                             */
/*------------------------------------------------------------------------*/
/*Description:  Provides control allocation functionality                 */
/*Type:         C - source file                                           */
/*Dependencies:                                                           */
/*------------------------------------------------------------------------*/
/*Author:       S. Hafner                                                 */
/*Date:         2023-11-08                                                */
/*Contact:      Institute of Flight System Dynamics                       */
/*              School of Engineering and Design                          */
/*              Technische Universitaet Muenchen                          */
/*              Boltzmannstrasse 15                                       */
/*              D-85748 Garching                                          */
/*              simon.hafner@tum.de                                       */
/*------------------------------------------------------------------------*/

#include <math.h>
#include "sym_mat_inv.h"
#include "control_allocation_LDL.h"

/**
 * @brief Does one control allocation step using the pseudo-inverse of B by LDL^T decomposition
 *
 * @param B Control effectiveness matrix
 * @param V Weighting vector with only 0 or positive entries
 * @param nu_des Desired pseudo control vector
 * @param u_alloc Allocated inputs u
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int CA_full_LDLT_six(matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, unsigned short *const rank, vector_t *u_alloc, const unsigned short rank_max)
{
    unsigned int status = 0;

    float A_matrix[6U * 6U] = {0.0f};
    matrix_t A =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = A_matrix};

    // TODO Do checks

    /* Scale B with weighting matrix */
    status |= mat_mult_diag_c(B, V);

    /* Calculate A */
    status |= mat_gram(B, &A);

    /* Run LDLT decomposition */
    status |= full_rank_LDLT(&A, L, d, rank, TOL, rank_max);

    // DO B^+ * nu
    status |= allocation_LDLT_six(B, V, nu_des, L, d, rank, u_alloc);

    return (status);
}

/**
 * @brief Simulink Interface Version:
 * Does one control allocation step using the pseudo-inverse of B by LDL^T decomposition
 *
 * @param num_u Number of inputs
 * @param B Control effectiveness matrix
 * @param V Weighting vector with only 0 or positive entries
 * @param nu_des Desired pseudo control vector
 * @param L Output L of the LDL^T
 * @param d Output d of the LDL^T
 * @param rank Rank of B
 * @param u_alloc Allocated inputs u
 * @param rank_max Maximum rank
 * @return unsigned int
 */
unsigned int sim_CA_full_LDLT_six(const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max)
{
    unsigned int status = 0;
    *rank = 6U;

    float A_matrix[6U * 6U] = {0.0f};
    matrix_t A =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = A_matrix};

    matrix_t B_m =
        {
            .rows = 6U,
            .columns = num_u,
            .matrix = B};

    vector_c_t V_M =
        {
            .rows = num_u,
            .vector = V};

    vector_c_t nu_des_m =
        {
            .rows = 6U,
            .vector = nu_des};

    matrix_t L_m =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L};

    vector_t d_m =
        {
            .rows = 6U,
            .vector = d};

    vector_t u_alloc_m =
        {
            .rows = num_u,
            .vector = u_alloc};

    // TODO Do checks

    mat_zero(&L_m);
    vec_zero(&d_m);

    /* Scale B with weighting matrix */
    status |= mat_mult_diag_c(&B_m, &V_M);

    /* Calculate A */
    status |= mat_gram(&B_m, &A);

    /* Run LDLT decomposition */
    status |= full_rank_LDLT(&A, &L_m, &d_m, rank, tol, rank_max);

    // DO B^+ * nu
    status |= allocation_LDLT_six(&B_m, &V_M, &nu_des_m, &L_m, &d_m, rank, &u_alloc_m);

    return (status);
}

/**
 * @brief Simulink Interface Version:
 * Does one control allocation step using the pseudo-inverse of B by LDL^T decomposition
 *
 * @param num_u Number of inputs
 * @param B Control effectiveness matrix
 * @param V Weighting vector with only 0 or positive entries
 * @param nu_des Desired pseudo control vector
 * @param L Output L of the LDL^T
 * @param d Output d of the LDL^T
 * @param rank Rank of B
 * @param u_alloc Allocated inputs u
 * @param rank_max Maximum rank
 * @return unsigned int
 */
unsigned int sim_CA_full_LDLT_four(const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max)
{
    unsigned int status = 0;
    *rank = 4U;

    float A_matrix[4U * 4U] = {0.0f};
    matrix_t A =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = A_matrix};

    matrix_t B_m =
        {
            .rows = 4U,
            .columns = num_u,
            .matrix = B};

    vector_c_t V_M =
        {
            .rows = num_u,
            .vector = V};

    vector_c_t nu_des_m =
        {
            .rows = 4U,
            .vector = nu_des};

    matrix_t L_m =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = L};

    vector_t d_m =
        {
            .rows = 4U,
            .vector = d};

    vector_t u_alloc_m =
        {
            .rows = num_u,
            .vector = u_alloc};

    // TODO Do checks

    mat_zero(&L_m);
    vec_zero(&d_m);

    /* Scale B with weighting matrix */
    status |= mat_mult_diag_c(&B_m, &V_M);

    /* Calculate A */
    status |= mat_gram(&B_m, &A);

    /* Run LDLT decomposition */
    status |= full_rank_LDLT(&A, &L_m, &d_m, rank, tol, rank_max);

    // DO B^+ * nu
    status |= allocation_LDLT_four(&B_m, &V_M, &nu_des_m, &L_m, &d_m, rank, &u_alloc_m);

    return (status);
}

/**
 * @brief Simulink Interface Version:
 * Does one control allocation step using the pseudo-inverse of B by pivoted LDL^T decomposition
 *
 * @param num_u Number of inputs
 * @param B Control effectiveness matrix
 * @param V Weighting vector with only 0 or positive entries
 * @param nu_des Desired pseudo control vector
 * @param L Output L of the LDL^T
 * @param d Output d of the LDL^T
 * @param rank Rank of B
 * @param u_alloc Allocated inputs u
 * @param rank_max Maximum rank
 * @return unsigned int
 */
unsigned int sim_CA_piv_LDLT_six(const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max)
{
    unsigned int status = 0;
    *rank = 6U;

    float A_matrix[6U * 6U] = {0.0f};
    unsigned short p_vector[6U] = {0U};
    matrix_t A =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = A_matrix};

    matrix_t B_m =
        {
            .rows = 6U,
            .columns = num_u,
            .matrix = B};

    vector_c_t V_M =
        {
            .rows = num_u,
            .vector = V};

    vector_c_t nu_des_m =
        {
            .rows = 6U,
            .vector = nu_des};

    matrix_t L_m =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L};

    vector_t d_m =
        {
            .rows = 6U,
            .vector = d};

    vector_t u_alloc_m =
        {
            .rows = num_u,
            .vector = u_alloc};

    vector_u_t p_m =
        {
            .rows = 6U,
            .vector = p_vector};

    // TODO Do checks

    mat_zero(&L_m);
    vec_zero(&d_m);

    /* Scale B with weighting matrix */
    status |= mat_mult_diag_c(&B_m, &V_M);

    /* Calculate A */
    status |= mat_gram(&B_m, &A);

    /* Run pivoted LDLT decomposition */
    status |= pivoted_LDLT(&A, &L_m, &d_m, &p_m, rank, tol, rank_max);

    // DO B^+ * nu
    status |= allocation_piv_LDLT_six(&B_m, &V_M, &nu_des_m, &L_m, &d_m, &p_m, rank, &u_alloc_m);

    return (status);
}

/**
 * @brief Simulink Interface Version:
 * Does one control allocation step using the pseudo-inverse of B by pivoted LDL^T decomposition
 *
 * @param num_u Number of inputs
 * @param B Control effectiveness matrix
 * @param V Weighting vector with only 0 or positive entries
 * @param nu_des Desired pseudo control vector
 * @param L Output L of the LDL^T
 * @param d Output d of the LDL^T
 * @param rank Rank of B
 * @param u_alloc Allocated inputs u
 * @param rank_max Maximum rank
 * @return unsigned int
 */
unsigned int sim_CA_piv_LDLT_four(const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max)
{
    unsigned int status = 0;
    *rank = 4U;

    float A_matrix[4U * 4U] = {0.0f};
    unsigned short p_vector[4U] = {0U};
    matrix_t A =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = A_matrix};

    matrix_t B_m =
        {
            .rows = 4U,
            .columns = num_u,
            .matrix = B};

    vector_c_t V_M =
        {
            .rows = num_u,
            .vector = V};

    vector_c_t nu_des_m =
        {
            .rows = 4U,
            .vector = nu_des};

    matrix_t L_m =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = L};

    vector_t d_m =
        {
            .rows = 4U,
            .vector = d};

    vector_t u_alloc_m =
        {
            .rows = num_u,
            .vector = u_alloc};

    vector_u_t p_m =
        {
            .rows = 4U,
            .vector = p_vector};

    // TODO Do checks

    mat_zero(&L_m);
    vec_zero(&d_m);

    /* Scale B with weighting matrix */
    status |= mat_mult_diag_c(&B_m, &V_M);

    /* Calculate A */
    status |= mat_gram(&B_m, &A);

    /* Run LDLT decomposition */
    status |= pivoted_LDLT(&A, &L_m, &d_m, &p_m, rank, tol, rank_max);

    // DO B^+ * nu
    status |= allocation_piv_LDLT_four(&B_m, &V_M, &nu_des_m, &L_m, &d_m, &p_m, rank, &u_alloc_m);

    return (status);
}

unsigned int allocation_LDLT_six(const matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, unsigned short *const rank, vector_t *u_alloc)
{
    unsigned int status = 0;

    float L_M_matrix[6U * 6U] = {0.0f};
    matrix_t L_M =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L_M_matrix};

    float M_matrix[6U * 6U] = {0.0f};
    matrix_t M =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = M_matrix};

    float L_sq_matrix[6U * 6U] = {0.0f};
    matrix_t L_sq =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L_sq_matrix};

    float L_M_inv_matrix[6U * 6U] = {0.0f};
    matrix_t L_M_inv =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L_M_inv_matrix};

    float d_M_vector[6U] = {0.0f};
    vector_t d_M = {6U, d_M_vector};
    float vec_inter_sc_vector[1U] = {0.0f};
    vector_t vec_inter_sc = {1U, vec_inter_sc_vector};
    float vec_inter_vector[6U] = {0.0f};
    vector_t vec_inter = {6U, vec_inter_vector};
    float vec_inter_2_vector[6U] = {0.0f};
    vector_t vec_inter_2 = {6U, vec_inter_2_vector};
    float vec_inter_3_vector[6U] = {0.0f};
    vector_t vec_inter_3 = {6U, vec_inter_3_vector};

    switch (*rank)
    {
    case 6U:
        /* B'*L_inv'*1/d*L_inv*nu_des */
        status |= solve_LDL_lower_vec_c(L, nu_des, &vec_inter);
        status |= diag_inv_vec(d, &vec_inter, &vec_inter_2);
        status |= solve_LDL_lower_T_vec(L, &vec_inter_2, &vec_inter_3);
        status |= mat_t_vec(B, &vec_inter_3, u_alloc);
        status |= diag_vec_c(V, u_alloc, u_alloc);

        break;

    case 5U:
        /* M inverse calculation */
        L->columns = *rank;
        L_sq.columns = *rank;
        L_sq.rows = *rank;
        L_M.columns = *rank;
        L_M.rows = *rank;
        L_M_inv.columns = *rank;
        L_M_inv.rows = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;
        d_M.rows = *rank;

        status |= lower_T_lower(L, &L_sq);

        /* Cholesky based inverse */
        // TODO Check if forward substitution is more efficient
        status |= full_rank_LDLT(&L_sq, &L_M, &d_M, rank, TOL, 5U);

        /* Invert diagonal & Lower*/
        status |= vec_inv(&d_M);
        status |= lower_inv_LDLT(&L_M, &L_M_inv);

        /* Calculate M */
        status |= mat_gram_T_weighted(&L_M_inv, &d_M, &M);

        if (status == 0U)
        {
            /* B'*L*M*1/d*M*L'*nu_des */
            vec_inter.rows = L->columns;
            status |= lower_t_vec_c(L, nu_des, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= lower_vec(L, &vec_inter_2, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }
        break;

    case 4U:
        /* M inverse calculation */
        L->columns = *rank;
        L_sq.columns = *rank;
        L_sq.rows = *rank;
        L_M.columns = *rank;
        L_M.rows = *rank;
        L_M_inv.columns = *rank;
        L_M_inv.rows = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;
        d_M.rows = *rank;

        status |= lower_T_lower(L, &L_sq);

        /* Cholesky based inverse */
        // TODO Check if forward substitution is more efficient
        status |= full_rank_LDLT(&L_sq, &L_M, &d_M, rank, TOL, 4U);

        /* Invert diagonal & Lower*/
        status |= vec_inv(&d_M);
        status |= lower_inv_LDLT(&L_M, &L_M_inv);

        /* Calculate M */
        status |= mat_gram_T_weighted(&L_M_inv, &d_M, &M);

        if (status == 0U)
        {
            /* B'*L*M*1/d*M*L'*nu_des */
            vec_inter.rows = L->columns;
            status |= lower_t_vec_c(L, nu_des, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= lower_vec(L, &vec_inter_2, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 3U:
        /* M inverse calculation */
        L->columns = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;

        status |= lower_T_lower(L, &L_sq);
        status |= sym_inv_three(&L_sq, &M);

        if (status == 0U)
        {
            /* B'*L*M*1/d*M*L'*nu_des */
            vec_inter.rows = L->columns;
            status |= lower_t_vec_c(L, nu_des, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= lower_vec(L, &vec_inter_2, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 2U:
        /* M inverse calculation */
        L->columns = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;

        status |= lower_T_lower(L, &L_sq);
        status |= sym_inv_two(&L_sq, &M);

        if (status == 0U)
        {
            /* B'*L*M*1/d*M*L'*nu_des */
            vec_inter.rows = L->columns;
            status |= lower_t_vec_c(L, nu_des, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= lower_vec(L, &vec_inter_2, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 1U:
        L->columns = *rank;

        status |= lower_T_lower(L, &L_sq);

        if ((ABS(L_sq.matrix[0]) > EPS_F) && (ABS(d->vector[0]) > EPS_F))
        {
            /* M inverse calculation */
            float M_sc = 1.0f / L_sq.matrix[0];

            /* B'*L*M*1/d*M*L'*nu_des */
            status |= lower_t_vec_c(L, nu_des, &vec_inter_sc);
            vec_scalar(&vec_inter_sc, (M_sc * (M_sc / d->vector[0])));
            status |= lower_vec(L, &vec_inter_sc, &vec_inter_2);
            status |= mat_t_vec(B, &vec_inter_2, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }
        else
        {
            status |= 2U;
        }
        break;

    case 0U:
        /* Return a zero vector for rank = 0 */
        break;

    default:
        /* Return a zero vector for rank = 0 */
        status |= 4U;
        break;
    }
    return (status);
}

unsigned int allocation_piv_LDLT_six(const matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, vector_u_t *p, unsigned short *const rank, vector_t *u_alloc)
{
    unsigned int status = 0;

    float L_M_matrix[6U * 6U] = {0.0f};
    matrix_t L_M =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L_M_matrix};

    float M_matrix[6U * 6U] = {0.0f};
    matrix_t M =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = M_matrix};

    float L_sq_matrix[6U * 6U] = {0.0f};
    matrix_t L_sq =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L_sq_matrix};

    float L_M_inv_matrix[6U * 6U] = {0.0f};
    matrix_t L_M_inv =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L_M_inv_matrix};

    float d_M_vector[6U] = {0.0f};
    vector_t d_M = {6U, d_M_vector};
    float vec_inter_sc_vector[1U] = {0.0f};
    vector_t vec_inter_sc = {1U, vec_inter_sc_vector};
    float vec_inter_vector[6U] = {0.0f};
    vector_t vec_inter = {6U, vec_inter_vector};
    float vec_inter_2_vector[6U] = {0.0f};
    vector_t vec_inter_2 = {6U, vec_inter_2_vector};
    float vec_inter_3_vector[6U] = {0.0f};
    vector_t vec_inter_3 = {6U, vec_inter_3_vector};

    switch (*rank)
    {
    case 6U:
        /* B'*P'*L_inv'*1/d*L_inv*P*nu_des */
        status |= vec_pivot_c(nu_des, p, &vec_inter);
        status |= solve_LDL_lower_vec(L, &vec_inter, &vec_inter_2);
        status |= diag_inv_vec(d, &vec_inter_2, &vec_inter);
        status |= solve_LDL_lower_T_vec(L, &vec_inter, &vec_inter_2);
        status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter_3);
        status |= mat_t_vec(B, &vec_inter_3, u_alloc);
        status |= diag_vec_c(V, u_alloc, u_alloc);

        break;

    case 5U:
        /* M inverse calculation */
        L->columns = *rank;
        L_sq.columns = *rank;
        L_sq.rows = *rank;
        L_M.columns = *rank;
        L_M.rows = *rank;
        L_M_inv.columns = *rank;
        L_M_inv.rows = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;
        d_M.rows = *rank;

        status |= lower_T_lower(L, &L_sq);

        /* Cholesky based inverse */
        // TODO Check if forward substitution is more efficient
        status |= full_rank_LDLT(&L_sq, &L_M, &d_M, rank, TOL, 5U);

        /* Invert diagonal & Lower*/
        status |= vec_inv(&d_M);
        status |= lower_inv_LDLT(&L_M, &L_M_inv);

        /* Calculate M */
        status |= mat_gram_T_weighted(&L_M_inv, &d_M, &M);

        if (status == 0U)
        {
            /* B'*P'*L*M*1/d*M*L'*P*nu_des */
            status |= vec_pivot_c(nu_des, p, &vec_inter);
            vec_inter_2.rows = L->columns;
            status |= lower_t_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = L->rows;
            status |= lower_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }
        break;

    case 4U:
        /* M inverse calculation */
        L->columns = *rank;
        L_sq.columns = *rank;
        L_sq.rows = *rank;
        L_M.columns = *rank;
        L_M.rows = *rank;
        L_M_inv.columns = *rank;
        L_M_inv.rows = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;
        d_M.rows = *rank;

        status |= lower_T_lower(L, &L_sq);

        /* Cholesky based inverse */
        // TODO Check if forward substitution is more efficient
        status |= full_rank_LDLT(&L_sq, &L_M, &d_M, rank, TOL, 4U);

        /* Invert diagonal & Lower*/
        status |= vec_inv(&d_M);
        status |= lower_inv_LDLT(&L_M, &L_M_inv);

        /* Calculate M */
        status |= mat_gram_T_weighted(&L_M_inv, &d_M, &M);

        if (status == 0U)
        {
            /* B'*P'*L*M*1/d*M*L'*P*nu_des */
            status |= vec_pivot_c(nu_des, p, &vec_inter);
            vec_inter_2.rows = L->columns;
            status |= lower_t_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = L->rows;
            status |= lower_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 3U:
        /* M inverse calculation */
        L->columns = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;

        status |= lower_T_lower(L, &L_sq);
        status |= sym_inv_three(&L_sq, &M);

        if (status == 0U)
        {
            /* B'*P'*L*M*1/d*M*L'*P*nu_des */
            status |= vec_pivot_c(nu_des, p, &vec_inter);
            vec_inter_2.rows = L->columns;
            status |= lower_t_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = L->rows;
            status |= lower_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 2U:
        /* M inverse calculation */
        L->columns = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;

        status |= lower_T_lower(L, &L_sq);
        status |= sym_inv_two(&L_sq, &M);

        if (status == 0U)
        {
            /* B'*P'*L*M*1/d*M*L'*P*nu_des */
            status |= vec_pivot_c(nu_des, p, &vec_inter);
            vec_inter_2.rows = L->columns;
            status |= lower_t_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = L->rows;
            status |= lower_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 1U:
        L->columns = *rank;

        status |= lower_T_lower(L, &L_sq);

        if ((ABS(L_sq.matrix[0]) > EPS_F) && (ABS(d->vector[0]) > EPS_F))
        {
            /* M inverse calculation */
            float M_sc = 1.0f / L_sq.matrix[0];

            /* B'*P'*L*M*1/d*M*L'*P*nu_des */
            status |= vec_pivot_c(nu_des, p, &vec_inter);
            status |= lower_t_vec(L, &vec_inter, &vec_inter_sc);
            vec_scalar(&vec_inter_sc, (M_sc * (M_sc / d->vector[0])));
            status |= lower_vec(L, &vec_inter_sc, &vec_inter_2);
            status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }
        else
        {
            status |= 2U;
        }
        break;

    case 0U:
        /* Return a zero vector for rank = 0 */
        break;

    default:
        /* Return a zero vector for rank = 0 */
        status |= 4U;
        break;
    }
    return (status);
}

unsigned int allocation_LDLT_four(const matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, unsigned short *const rank, vector_t *u_alloc)
{
    unsigned int status = 0;

    float M_matrix[4U * 4U] = {0.0f};
    matrix_t M =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = M_matrix};

    float L_sq_matrix[4U * 4U] = {0.0f};
    matrix_t L_sq =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = L_sq_matrix};

    float vec_inter_sc_vector[1U] = {0.0f};
    vector_t vec_inter_sc = {1U, vec_inter_sc_vector};
    float vec_inter_vector[4U] = {0.0f};
    vector_t vec_inter = {4U, vec_inter_vector};
    float vec_inter_2_vector[4U] = {0.0f};
    vector_t vec_inter_2 = {4U, vec_inter_2_vector};
    float vec_inter_3_vector[4U] = {0.0f};
    vector_t vec_inter_3 = {4U, vec_inter_3_vector};

    switch (*rank)
    {
    case 4U:
        /* B'*L_inv'*1/d*L_inv*nu_des */
        status |= solve_LDL_lower_vec_c(L, nu_des, &vec_inter);
        status |= diag_inv_vec(d, &vec_inter, &vec_inter_2);
        status |= solve_LDL_lower_T_vec(L, &vec_inter_2, &vec_inter_3);
        status |= mat_t_vec(B, &vec_inter_3, u_alloc);
        status |= diag_vec_c(V, u_alloc, u_alloc);

        break;

    case 3U:
        /* M inverse calculation */
        L->columns = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;

        status |= lower_T_lower(L, &L_sq);
        status |= sym_inv_three(&L_sq, &M);

        if (status == 0U)
        {
            /* B'*L*M*1/d*M*L'*nu_des */
            vec_inter.rows = L->columns;
            status |= lower_t_vec_c(L, nu_des, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= lower_vec(L, &vec_inter_2, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 2U:
        /* M inverse calculation */
        L->columns = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;

        status |= lower_T_lower(L, &L_sq);
        status |= sym_inv_two(&L_sq, &M);

        if (status == 0U)
        {
            /* B'*L*M*1/d*M*L'*nu_des */
            vec_inter.rows = L->columns;
            status |= lower_t_vec_c(L, nu_des, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= mat_vec(&M, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= lower_vec(L, &vec_inter_2, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 1U:
        L->columns = *rank;

        status |= lower_T_lower(L, &L_sq);

        if ((ABS(L_sq.matrix[0]) > EPS_F) && (ABS(d->vector[0]) > EPS_F))
        {
            /* M inverse calculation */
            float M_sc = 1.0f / L_sq.matrix[0];

            /* B'*L*M*1/d*M*L'*nu_des */
            status |= lower_t_vec_c(L, nu_des, &vec_inter_sc);
            vec_scalar(&vec_inter_sc, (M_sc * (M_sc / d->vector[0])));
            status |= lower_vec(L, &vec_inter_sc, &vec_inter_2);
            status |= mat_t_vec(B, &vec_inter_2, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }
        else
        {
            status |= 2U;
        }
        break;

    case 0U:
        /* Return a zero vector for rank = 0 */
        break;

    default:
        /* Return a zero vector for rank = 0 */
        status |= 4U;
        break;
    }
    return (status);
}

unsigned int allocation_piv_LDLT_four(const matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, vector_u_t *p, unsigned short *const rank, vector_t *u_alloc)
{
    unsigned int status = 0;

    float M_matrix[4U * 4U] = {0.0f};
    matrix_t M =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = M_matrix};

    float L_sq_matrix[4U * 4U] = {0.0f};
    matrix_t L_sq =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = L_sq_matrix};

    float vec_inter_sc_vector[1U] = {0.0f};
    vector_t vec_inter_sc = {1U, vec_inter_sc_vector};
    float vec_inter_vector[4U] = {0.0f};
    vector_t vec_inter = {4U, vec_inter_vector};
    float vec_inter_2_vector[4U] = {0.0f};
    vector_t vec_inter_2 = {4U, vec_inter_2_vector};
    float vec_inter_3_vector[4U] = {0.0f};
    vector_t vec_inter_3 = {4U, vec_inter_3_vector};

    switch (*rank)
    {
    case 4U:
        /* B'*P'*L_inv'*1/d*L_inv*P*nu_des */
        status |= vec_pivot_c(nu_des, p, &vec_inter);
        status |= solve_LDL_lower_vec(L, &vec_inter, &vec_inter_2);
        status |= diag_inv_vec(d, &vec_inter_2, &vec_inter);
        status |= solve_LDL_lower_T_vec(L, &vec_inter, &vec_inter_2);
        status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter_3);
        status |= mat_t_vec(B, &vec_inter_3, u_alloc);
        status |= diag_vec_c(V, u_alloc, u_alloc);

        break;

    case 3U:
        /* M inverse calculation */
        L->columns = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;

        status |= lower_T_lower(L, &L_sq);
        status |= sym_inv_three(&L_sq, &M);

        if (status == 0U)
        {
            /* B'*P'*L*M*1/d*M*L'*P*nu_des */
            status |= vec_pivot_c(nu_des, p, &vec_inter);
            vec_inter_2.rows = L->columns;
            status |= lower_t_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = L->rows;
            status |= lower_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 2U:
        /* M inverse calculation */
        L->columns = *rank;
        M.columns = *rank;
        M.rows = *rank;
        d->rows = *rank;

        status |= lower_T_lower(L, &L_sq);
        status |= sym_inv_two(&L_sq, &M);

        if (status == 0U)
        {
            /* B'*P'*L*M*1/d*M*L'*P*nu_des */
            status |= vec_pivot_c(nu_des, p, &vec_inter);
            vec_inter_2.rows = L->columns;
            status |= lower_t_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = M.columns;
            status |= diag_inv_vec(d, &vec_inter, &vec_inter_2);
            vec_inter.rows = M.columns;
            status |= mat_vec(&M, &vec_inter_2, &vec_inter);
            vec_inter_2.rows = L->rows;
            status |= lower_vec(L, &vec_inter, &vec_inter_2);
            vec_inter.rows = L->rows;
            status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }

        break;

    case 1U:
        L->columns = *rank;

        status |= lower_T_lower(L, &L_sq);

        if ((ABS(L_sq.matrix[0]) > EPS_F) && (ABS(d->vector[0]) > EPS_F))
        {
            /* M inverse calculation */
            float M_sc = 1.0f / L_sq.matrix[0];

            /* B'*P'*L*M*1/d*M*L'*P*nu_des */
            status |= vec_pivot_c(nu_des, p, &vec_inter);
            status |= lower_t_vec(L, &vec_inter, &vec_inter_sc);
            vec_scalar(&vec_inter_sc, (M_sc * (M_sc / d->vector[0])));
            status |= lower_vec(L, &vec_inter_sc, &vec_inter_2);
            status |= vec_pivot_rev(&vec_inter_2, p, &vec_inter);
            status |= mat_t_vec(B, &vec_inter, u_alloc);
            status |= diag_vec_c(V, u_alloc, u_alloc);
        }
        else
        {
            status |= 2U;
        }
        break;

    case 0U:
        /* Return a zero vector for rank = 0 */
        break;

    default:
        /* Return a zero vector for rank = 0 */
        status |= 4U;
        break;
    }
    return (status);
}

unsigned int CA_down_LDLT_six(const unsigned short saturated, matrix_t *const B, const vector_c_t *const V, const vector_c_t *const nu_des, matrix_t *L, vector_t *d, unsigned short *const rank, vector_t *u_alloc, const unsigned short rank_max)
{
    unsigned int status = 0;
    float v_vector[6U] = {0.0f};
    vector_t v = {6U, v_vector};

    /* Select the correct column of B as v*/
    if (saturated >= B->columns)
    {
        status |= 1U;
    }
    else
    {
        for (unsigned short i = 0; i < B->rows; i++)
        {
            /* Need square root of V as we defined A = BWB^T*/
            v.vector[i] = M_P_ENTRY(B, i, saturated);
        }
    }

    /* Set B of saturated to 0 */
    for (unsigned short i = 0; i < B->rows; i++)
    {
        /* Need square root of V as we defined A = BWB^T*/
        M_P_ENTRY(B, i, saturated) = 0.0f;
    }

    status |= rank_one_LDLT_down(v, L, d, rank, 12.0f * TOL, rank_max);

    status |= allocation_LDLT_six(B, V, nu_des, L, d, rank, u_alloc);

    /* Square Weighting Vector */
    status |= diag_vec_c(V, u_alloc, u_alloc);

    return (status);
}

/**
 * @brief Simulink interface version of the staturation rank one step
 *
 * @param saturated Saturated input of this step
 * @param num_u Number of inputs
 * @param B Control effectiveness matrix
 * @param V Weighting vector V = V*V
 * @param nu_des Desired pseudo control vector
 * @param L L of LDL^T to be updated
 * @param d d of LDL^T to be updated
 * @param rank Rank of B also input
 * @param u_alloc Allocated inputs u
 * @param rank_max Maximum rank
 * @return unsigned int
 */
unsigned int sim_CA_down_LDLT_six(const unsigned short saturated, const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max)
{
    unsigned int status = 0;
    float v_vector[6U] = {0.0f};
    vector_t v = {6U, v_vector};

    matrix_t B_m =
        {
            .rows = 6U,
            .columns = num_u,
            .matrix = B};

    vector_c_t V_m =
        {
            .rows = num_u,
            .vector = V};

    vector_c_t nu_des_m =
        {
            .rows = 6U,
            .vector = nu_des};

    matrix_t L_m =
        {
            .rows = 6U,
            .columns = 6U,
            .matrix = L};

    vector_t d_m =
        {
            .rows = 6U,
            .vector = d};

    vector_t u_alloc_m =
        {
            .rows = num_u,
            .vector = u_alloc};

    /* Select the correct column of B as v*/
    if (saturated >= B_m.columns)
    {
        status |= 1U;
    }
    else
    {
        for (unsigned short i = 0; i < B_m.rows; i++)
        {
            /* Extract Vector from B*/
            v.vector[i] = M_ENTRY(B_m, i, saturated);
        }
    }

    /* Set B of saturated to 0 */
    for (unsigned short i = 0; i < B_m.rows; i++)
    {
        M_ENTRY(B_m, i, saturated) = 0.0f;
    }

    status |= rank_one_LDLT_down(v, &L_m, &d_m, rank, tol, rank_max);

    status |= allocation_LDLT_six(&B_m, &V_m, &nu_des_m, &L_m, &d_m, rank, &u_alloc_m);

    /* Ensure saturated u_alloc_element is exactly 0*/
    u_alloc_m.vector[saturated] = 0.0f;

    return (status);
}

/**
 * @brief Simulink interface version of the staturation rank one step
 *
 * @param saturated Saturated input of this step
 * @param num_u Number of inputs
 * @param B Control effectiveness matrix
 * @param V Weighting vector V = V*V
 * @param nu_des Desired pseudo control vector
 * @param L L of LDL^T to be updated
 * @param d d of LDL^T to be updated
 * @param rank Rank of B also input
 * @param u_alloc Allocated inputs u
 * @param rank_max Maximum rank
 * @return unsigned int
 */
unsigned int sim_CA_down_LDLT_four(const unsigned short saturated, const unsigned short num_u, const float tol, float *const B, const float *const V, const float *const nu_des, float *const L, float *const d, unsigned short *const rank, float *const u_alloc, const unsigned short rank_max)
{
    unsigned int status = 0;
    float v_vector[4U] = {0.0f};
    vector_t v = {4U, v_vector};

    matrix_t B_m =
        {
            .rows = 4U,
            .columns = num_u,
            .matrix = B};

    vector_c_t V_m =
        {
            .rows = num_u,
            .vector = V};

    vector_c_t nu_des_m =
        {
            .rows = 4U,
            .vector = nu_des};

    matrix_t L_m =
        {
            .rows = 4U,
            .columns = 4U,
            .matrix = L};

    vector_t d_m =
        {
            .rows = 4U,
            .vector = d};

    vector_t u_alloc_m =
        {
            .rows = num_u,
            .vector = u_alloc};

    /* Select the correct column of B as v*/
    if (saturated >= B_m.columns)
    {
        status |= 1U;
    }
    else
    {
        for (unsigned short i = 0; i < B_m.rows; i++)
        {
            /* Need square root of V equal to V as we defined A = BWB^T*/
            v.vector[i] = M_ENTRY(B_m, i, saturated);
        }
    }

    /* Set B of saturated to 0 */
    for (unsigned short i = 0; i < B_m.rows; i++)
    {
        /* Need square root of V as we defined A = BWB^T*/
        M_ENTRY(B_m, i, saturated) = 0.0f;
    }

    status |= rank_one_LDLT_down(v, &L_m, &d_m, rank, tol, rank_max);

    status |= allocation_LDLT_four(&B_m, &V_m, &nu_des_m, &L_m, &d_m, rank, &u_alloc_m);

    /* Ensure saturated u_alloc_element is exactly 0*/
    u_alloc_m.vector[saturated] = 0.0f;

    return (status);
}

/**
 * @brief Calculates the full rank Cholesky decomposition $A = LL^T$ based on
 * P. Courrieu, "Fast Computation of Moore-Penrose Inverse Matrices", Neural Information Processing, 2005.
 * arxiv.org/abs/0804.4809
 *
 * @param A Squared positive semidefinite Matrix
 * @param L Lower matrix result of the Cholesky decomposition
 * @param rank rank of the matrix A and number of non-zero columns in L
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int full_rank_Cholesky(const matrix_t *const A, matrix_t *L, unsigned short *const rank)
{
    unsigned int status = 0;
    if (A->rows != A->columns)
    {
        status = 1U; /* return status 1 (ERROR) Matrix non - squared */
    }
    else
    {
        float max_diag = mat_max_diag(A);
        float tol = (float)A->rows * EPS_F * max_diag;

        unsigned short r = 0U;

        for (unsigned short k = 0U; k < A->rows; k++)
        {
            if (r == 0U)
            {
                for (unsigned short i = k; i < A->rows; i++)
                {
                    M_P_ENTRY(L, i, r) = M_P_ENTRY(A, i, k);
                }
            }
            else
            {
                for (unsigned short i = k; i < A->rows; i++)
                {
                    // Compute inner product LL = L[k:n, :r] * L[k, :r].T
                    float LL = 0.0f;

                    for (unsigned short j = 0U; j < r; j++)
                    {
                        LL += M_P_ENTRY(L, i, j) * M_P_ENTRY(L, k, j);
                    }

                    M_P_ENTRY(L, i, r) = M_P_ENTRY(A, i, k) - LL;
                }
            }

            if (M_P_ENTRY(L, k, r) > tol)
            {
                M_P_ENTRY(L, k, r) = sqrtf(M_P_ENTRY(L, k, r)); // TODO Use square root of Kajetan Nürnberger

                if (k < (A->rows - 1U))
                {
                    for (unsigned short i = k + 1U; i < A->rows; i++)
                    {
                        M_P_ENTRY(L, i, r) = M_P_ENTRY(L, i, r) / M_P_ENTRY(L, k, r);
                    }
                }

                r += 1U;
            }
        }
        *rank = r; /* Return Rank */
    } /* return status 0 (OK) */
    return (status);
}

/**
 * @brief Conducts the full-rank LDL^T decomposition based on the cholesky decomposition by
 * P. Courrieu, “Fast Computation of Moore-Penrose Inverse Matrices,” Neural Information Processing, 2005.
 *
 * @param A Symmetric, positive semidefinite matrix (n x n), with rank rank
 * @param L Resulting Lower matrix with ones on the diagonal
 * @param d Vector containing the diagonal values of D
 * @param rank Rank of A
 * @param rank_max Theoretically maximum rank of A
 * @return unsigned int 0 = successful exit; 1 = wrong dimensions
 */
unsigned int full_rank_LDLT(const matrix_t *const A, matrix_t *L, vector_t *d, unsigned short *const rank, const float tol, const unsigned short rank_max)
{
    unsigned int status = 0;

    if ((A->rows != A->columns) || (A->rows != d->rows))
    {
        status = 1U; /* return status 1 (ERROR) Matrix non - squared */
    }
    else
    {
        float max_diag = mat_max_diag(A);
        float tol_A = tol * ABS(max_diag);

        unsigned short r = 0U;

        // TODO store D within L
        for (unsigned short k = 0U; k < A->rows; k++)
        {
            if (r == 0U)
            {
                d->vector[r] = M_P_ENTRY(A, k, k);
                for (unsigned short i = (k + 1U); i < A->rows; i++)
                {
                    M_P_ENTRY(L, i, r) = M_P_ENTRY(A, i, k);
                }
            }
            else
            {
                // Compute diagonal element
                float LL = 0.0f;

                for (unsigned short j = 0U; j < r; j++)
                {
                    LL += M_P_ENTRY(L, k, j) * (M_P_ENTRY(L, k, j) * d->vector[j]);
                }
                d->vector[r] = M_P_ENTRY(A, k, k) - LL;

                for (unsigned short i = (k + 1U); i < A->rows; i++)
                {
                    // A((k+1):n,k) - L((k+1):n, 1:(r-1))*(L(k,1:(r-1))'.*d(1:(r-1)))
                    LL = 0.0f;

                    for (unsigned short j = 0U; j < r; j++)
                    {
                        LL += M_P_ENTRY(L, i, j) * (M_P_ENTRY(L, k, j) * d->vector[j]);
                    }

                    M_P_ENTRY(L, i, r) = M_P_ENTRY(A, i, k) - LL;
                }
            }
            if ((d->vector[r] > tol_A) && (r < rank_max))
            {
                M_P_ENTRY(L, k, r) = 1.0f;
                for (unsigned short i = (k + 1U); i < A->rows; i++)
                {
                    M_P_ENTRY(L, i, r) = M_P_ENTRY(L, i, r) / d->vector[r];
                }

                r += 1U;
            }
        }
        *rank = r; /* Return Rank */
    } /* return status 0 (OK) */
    return (status);
}

unsigned int pivoted_LDLT(const matrix_t *const A, matrix_t *L, vector_t *d, vector_u_t *p, unsigned short *const rank, const float tol, const unsigned short rank_max)
{
    unsigned int status = 0;

    if ((A->rows != A->columns) || (A->rows != d->rows))
    {
        status = 1U; /* return status 1 (ERROR) Matrix non - squared */
    }
    else
    {
        float max_diag = mat_max_diag(A);
        float tol_A = tol * ABS(max_diag);

        unsigned short r = 0U;

        /*Initialize p*/
        for (unsigned short i = 0; i < p->rows; i++)
        {
            p->vector[i] = i;
        }

        for (unsigned short k = 0U; k < A->rows; k++)
        {
            unsigned short idx_pivot = mat_max_diag_pivot(A, p, k);
            unsigned short piv_k = p->vector[k];
            p->vector[k] = p->vector[idx_pivot];
            p->vector[idx_pivot] = piv_k;

            unsigned short p_k = p->vector[k];

            /*Assigning d out of pivoted sub A*/
            d->vector[k] = M_P_ENTRY(A, p_k, p_k);

            if ((d->vector[k] > tol_A) && (r < rank_max))
            {
                /*A(piv((k+1):n),piv((k+1):n)) = A(piv((k+1):n),piv((k+1):n)) - (A(piv((k+1):n),piv(k))*A(piv((k+1):n),piv(k))')/d(k)*/
                for (unsigned short i = (k + 1U); i < A->rows; i++)
                {
                    unsigned short p_i = p->vector[i];
                    for (unsigned short j = (k + 1U); j < A->rows; j++)
                    {
                        M_P_ENTRY(A, p_i, p->vector[j]) -= (M_P_ENTRY(A, p_i, p_k) * M_P_ENTRY(A, p->vector[j], p_k)) / d->vector[k];
                    }
                }
                /*A(piv((k+1):n),piv(k)) = A(piv((k+1):n),piv(k))/d(k)*/
                for (unsigned short i = (k + 1U); i < A->rows; i++)
                {
                    M_P_ENTRY(A, p->vector[i], p_k) /= d->vector[k];
                }
                r += 1U;
            }
            else
            {
                d->vector[k] = 0.0f;
            }
        }

        /*Copy Resulting L*/
        // TODO Can be avoided if we directly use A
        for (unsigned short k = 0; k < r; k++)
        {
            unsigned short p_k = p->vector[k];
            M_P_ENTRY(L, k, k) = 1.0f;
            for (unsigned short i = (k + 1U); i < A->rows; i++)
            {
                M_P_ENTRY(L, i, k) = M_P_ENTRY(A, p->vector[i], p_k);
            }
        }

        *rank = r; /* Return Rank */
    } /* return status 0 (OK) */
    return (status);
}

/**
 * @brief Performs a rank one down date of the matrix $A = L*D*L^T$ by the vector v such that $A' = A - c*c^T$
 *
 * @param v Rank one update vector
 * @param L Lower matrix of the LDL'T of A. Is overwritten
 * @param d Vector with the diagonal elements of D
 * @param rank Current rank of L
 * @return unsigned int 0 = OK
 */
unsigned int rank_one_LDLT_down(vector_t v, matrix_t *L, vector_t *d, unsigned short *const rank, const float tol, const unsigned short rank_max)
{
    unsigned int status = 0U;
    unsigned short matr_row_size = L->rows;
    float p = 0.0f;
    float alpha = -1.0f;
    unsigned short r = 0U;

    if (*rank > matr_row_size)
    {
        status = 1U;
    }

    float max_diag = vec_max(d);
    float tol_A = tol * ABS(max_diag);

    for (unsigned short j = 0U; j < *rank; j++)
    {
        for (unsigned short i = j; i < matr_row_size; i++)
        {
            if ((M_P_ENTRY(L, i, j) < (1.0f + EPS_F)) && (M_P_ENTRY(L, i, j) > (1.0f - EPS_F)))
            {
                p = v.vector[i];
                break;
            }
        }
        for (unsigned short i = (j + 1U); i < matr_row_size; i++)
        {
            v.vector[i] = v.vector[i] - (p * M_P_ENTRY(L, i, j));
        }

        float d_j = d->vector[j];
        d->vector[r] = d_j + (alpha * (p * p));
        if ((d->vector[r] > tol_A) && (r < rank_max))
        {
            float beta = (p * alpha) / d->vector[r];
            alpha = (d_j * alpha) / d->vector[r];
            for (unsigned short i = (j + 1U); i < matr_row_size; i++)
            {
                M_P_ENTRY(L, i, r) = M_P_ENTRY(L, i, j) + (beta * v.vector[i]);
            }
            M_P_ENTRY(L, j, r) = 1.0f;
            r += 1U;
        }
        else
        {
            d->vector[r] = 0.0f;
            for (unsigned short i = j; i < matr_row_size; i++)
            {
                M_P_ENTRY(L, i, r) = 0.0f;
            }
        }
    }
    *rank = r;
    return (status);
}