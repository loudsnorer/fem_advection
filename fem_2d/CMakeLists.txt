xilinx_synthesis_target("entry" "" "solve_jacobi" "jacobi_solver_stencil.cpp")

xilinx_testbench("fem_stencil.cpp jacobi_solver_stencil entry.cpp" "entry")
