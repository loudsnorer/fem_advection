#include "jacobi_solver_stencil.h"

void entry(double stencil[stencil_size], int indirections[stencil_size*nodes_number], double b[nodes_number], double solution_1[nodes_number], double solution_2[nodes_number])
{
    /*
    #pragma HLS INTERFACE m_axi port=indirections bundle=gmem0 offset=slave
    #pragma HLS INTERFACE m_axi port=b bundle=gmem0 offset=slave
    #pragma HLS INTERFACE m_axi port=solution1 bundle=gmem0 offset=slave
    #pragma HLS INTERFACE m_axi port=solution2 bundle=gmem1 offset=slave
    #pragma HLS INTERFACE s_axilite port=indirections
    #pragma HLS INTERFACE s_axilite port=b
    #pragma HLS INTERFACE s_axilite port=solution1
    #pragma HLS INTERFACE s_axilite port=solution2
    #pragma HLS INTERFACE s_axilite port=return
    */
    solve_jacobi(stencil, indirections, b, solution_1, solution_2);

}