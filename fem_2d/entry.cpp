#include "jacobi_solver_stencil.h"

void entry(int num_iters, double stencil[9], int* indirections, double* b, int nodes_number, int nodes_number_x, int nodes_number_y, double* solution_1, double* solution_2)
{
    #pragma HLS INTERFACE m_axi port=indirections bundle=gmem0 offset=slave
    #pragma HLS INTERFACE m_axi port=b bundle=gmem0 offset=slave
    #pragma HLS INTERFACE m_axi port=solution1 bundle=gmem0 offset=slave
    #pragma HLS INTERFACE m_axi port=solution2 bundle=gmem0 offset=slave
    #pragma HLS INTERFACE s_axilite port=indirections
    #pragma HLS INTERFACE s_axilite port=b
    #pragma HLS INTERFACE s_axilite port=solution1
    #pragma HLS INTERFACE s_axilite port=solution2
    #pragma HLS INTERFACE s_axilite port=return

    solve_jacobi(num_iters, stencil, indirections, b,nodes_number, nodes_number_x, nodes_number_y, solution_1, solution_2);

}