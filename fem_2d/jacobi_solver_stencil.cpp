#include "jacobi_solver_stencil.h"
#include <iostream>

void solve_jacobi(double stencil[stencil_size], int indirections[stencil_size*nodes_number], double b[nodes_number], double solution_1[nodes_number], double solution_2[nodes_number])
{    
    const auto st0 = stencil[0];
    const auto st1 = stencil[1];
    const auto st2 = stencil[2];
    const auto st3 = stencil[3];
    const auto st4 = stencil[4];
    const auto st5 = stencil[5];
    const auto st6 = stencil[6];
    const auto st7 = stencil[7];
    const auto st8 = stencil[8];

    for (int i = 0; i < nodes_number; i++)
    {
        solution_1[i] = 0.0;
    }

    // iterations of jacobi
    // simple jacobi solver
    // note Dirichlet BC is implied
    double *s1 = solution_1;
    double *s2 = solution_2;
    double *tmp;
    for (size_t it = 0; it < num_iters; it++)
    {

         // The first two rows are buffered in separate pipelines

        double V0[nodes_number_x];
        double V1[nodes_number_x];

        for (size_t j = 0; j < nodes_number_x; j++)
        {
            #pragma HLS pipeline
            V0[j] = s1[j];
        }

        for (size_t j = 0; j < nodes_number_x; j++)
        {
            #pragma HLS pipeline
            V1[j] = s1[j+nodes_number_x];
        }

        // The remaining rows can be streamed


        for (size_t idx_i = 1; idx_i < nodes_number_x - 1; idx_i++)
        {
            for (size_t idx_j = 1; idx_j < nodes_number_y - 1; idx_j++)
            {
				#pragma HLS pipeline II=1
                #pragma HLS loop_flatten
                int i = idx_i * nodes_number_y + idx_j;

                /*
                const auto v0 = s1[indirections[stencil_size * i + 0]];
                const auto v1 = s1[indirections[stencil_size * i + 1]];
                const auto v2 = s1[indirections[stencil_size * i + 2]];
                const auto v3 = s1[indirections[stencil_size * i + 3]];
                const auto v4 = s1[indirections[stencil_size * i + 4]];
                const auto v5 = s1[indirections[stencil_size * i + 5]];
                const auto v6 = s1[indirections[stencil_size * i + 6]];
                const auto v7 = s1[indirections[stencil_size * i + 7]];
                const auto v8 = s1[indirections[stencil_size * i + 8]];
                const auto rhs = b[i];
                */

                const auto V2_0 = s1[indirections[stencil_size * i + 6]];
                const auto V2_1 = s1[indirections[stencil_size * i + 7]];
                const auto V2_2 = s1[indirections[stencil_size * i + 8]];
                const auto rhs = b[i];

                const auto result = 1.0 / st4 *
                             (rhs - 
                             (st0 * V0[idx_j-1] + 
                             st1 * V0[idx_j  ] + 
                             st2 * V0[idx_j+1] + 
                             st3 * V1[idx_j-1] +
                             st5 * V1[idx_j+1] + 
                             st6 * V2_0 + 
                             st7 * V2_1 + 
                             st8 * V2_2));

                // we are ignoring BC here, but it should be considered in general case
                V0[idx_j]=V1[idx_j];
                V1[idx_j]=V2_1;

                #pragma HLS dependence variable=V0 false
                #pragma HLS dependence variable=V1 false

                s2[i] = result;
            }
        }

        // swap the buffers
        tmp = s1;
        s1 = s2;
        s2 = tmp;
    }
}
