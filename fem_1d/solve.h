//
// Created by Fariz Huseynli on 12.11.20.
//

#ifndef FEM_ADVECTION_SOLVE_H
#define FEM_ADVECTION_SOLVE_H

#endif //FEM_ADVECTION_SOLVE_H

double **r8rmat_copy_new ( int m, int n, double **a );
void r8rmat_delete ( int m, int n, double **a );
double *r8rmat_fs_new ( int n, double **a, double b[] );
double **r8rmat_new ( int m, int n );
double **r8rmat_zero ( int m, int n );
double *r8vec_copy_new ( int n, double a1[] );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );
