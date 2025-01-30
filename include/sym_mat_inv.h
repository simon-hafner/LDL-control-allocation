/*------------------------------------------------------------------------*/
/*                           sym_mat_inv.h                                */
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
/*Description:  Provides inverses of symmetric matrices size 2 to 6       */
/*Type:         C - header file                                           */
/*Dependencies:                                                           */
/*------------------------------------------------------------------------*/
/*Author:       S. Hafner                                                 */
/*Date:         2023-10-25                                                */
/*Contact:      Institute of Flight System Dynamics                       */
/*              School of Engineering and Design                          */
/*              Technische Universitaet Muenchen                          */
/*              Boltzmannstrasse 15                                       */
/*              D-85748 Garching                                          */
/*              simon.hafner@tum.de                                       */
/*------------------------------------------------------------------------*/

#ifndef SYM_MAT_INV_H
#define SYM_MAT_INV_H

#include "mat_vec.h"

unsigned int sym_inv_two(const matrix_t *const matr, matrix_t *res_matr);

unsigned int sym_inv_three(const matrix_t *const matr, matrix_t *res_matr);

unsigned int sym_inv_four(const matrix_t *const matr, matrix_t *res_matr);

unsigned int sym_inv_five(const matrix_t *const matr, matrix_t *res_matr);

// unsigned int sym_inv_six(matrix_t matr, matrix_t *res_matr);

#endif /*SYM_MAT_INV_H*/