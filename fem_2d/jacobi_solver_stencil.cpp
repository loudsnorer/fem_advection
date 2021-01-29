#include "jacobi_solver_stencil.h"
#include <iostream>

void solve_jacobi(double stencil[stencil_size], int indirections[stencil_size*nodes_number], double b[nodes_number], double solution_1[nodes_number], double solution_2[nodes_number])
{    
    for (int i = 0; i < nodes_number; i++)
    {
        solution_1[i] = 0.0;
    }

    // iterations of jacobi
    double *s1 = solution_1;
    double *s2 = solution_2;
    double *tmp;
    for (size_t it = 0; it < num_iters; it++)
    {
        // simple jacobi solver
        // note Dirichlet BC is implied
        for (size_t idx_i = 1; idx_i < nodes_number_x - 1; idx_i++)
        {
            for (size_t idx_j = 1; idx_j < nodes_number_y - 1; idx_j++)
            {
				#pragma HLS pipeline II=1
                int i = idx_i * nodes_number_y + idx_j;
                s2[i] = 1.0 / stencil[4] *
                             (b[i] - 
                             (stencil[0] * s1[indirections[stencil_size * i + 0]] + 
                             stencil[1] * s1[indirections[stencil_size * i + 1]] + 
                             stencil[2] * s1[indirections[stencil_size * i + 2]] + 
                             stencil[3] * s1[indirections[stencil_size * i + 3]] + 
                             stencil[5] * s1[indirections[stencil_size * i + 5]] + 
                             stencil[6] * s1[indirections[stencil_size * i + 6]] + 
                             stencil[7] * s1[indirections[stencil_size * i + 7]] + 
                             stencil[8] * s1[indirections[stencil_size * i + 8]]));
            }
        }

        // swap the buffers
        tmp = s1;
        s1 = s2;
        s2 = tmp;
    }
}
