#pragma once

#define PI 3.14159265358979323846
#define nodes_number_x 16
#define nodes_number_y 16
#define nodes_number nodes_number_x * nodes_number_y
#define num_iters 1
#define stencil_size 9

void solve_jacobi(double stencil[stencil_size], int indirections[stencil_size*nodes_number], double b[nodes_number], double solution_1[nodes_number], double solution_2[nodes_number]);
void entry(double stencil[stencil_size], int indirections[stencil_size*nodes_number], double b[nodes_number], double solution_1[nodes_number], double solution_2[nodes_number]);