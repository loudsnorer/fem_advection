#include <iostream>
#include <math.h>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
#include <algorithm>
#include <chrono>
#include <iostream>
#include<vector>
#include <fstream>

using namespace std;
using namespace std::chrono;

using namespace std;

# include "solve.h"
#define PI 3.14159265
constexpr int nodes_number{ 4096    };

double fRand(double i, int i1);
void create_matrix();
double calculate_fi(double x, double *nodes, int i1);

double a_sin_bx(double param1, double param2, double x);
void create_nodes (double* nodes);
void create_matrix_b(double *b, double *f);
double calculate_dfi(double x, double *nodes, int i);
void create_matrix_K();
double * solve_matrix();
void check_solution(double *pDouble);
void write_solutions_to_file();



ofstream command_unit;
string data_filename = "fem_data.3d";
ofstream data_unit;



//number of points to be calculated for each function
//TO-DO: implement in a decent way
constexpr int iter{ 1 };


//pointer to main function values calculated for matrix
double f_x[nodes_number][iter]{};

//pointer to main function values calculated for matrix in 1D
double f_x_1d[nodes_number + 2];

//pointer to values of base functions calculated for matrix in 1D
double **fi_p;
double fi[nodes_number + 1][nodes_number + 2];
double dfi[nodes_number + 1][nodes_number + 2];

double nodes[nodes_number + 2];

//b vector
double b[nodes_number];

//K matrix
double K[nodes_number][nodes_number]{};

double* solution;

int main() {
//    std::cout << "Hello, World!" << std::endl;

    // Get starting timepoint
    auto start_generate = high_resolution_clock::now();
    create_matrix();


    create_matrix_b(b, f_x_1d);

    create_matrix_K();

    // Get ending timepoint
    auto stop_generate = high_resolution_clock::now();




    // Get starting timepoint
    auto start_solve = high_resolution_clock::now();
    solution = solve_matrix();

    // Get ending timepoint
    auto stop_solve = high_resolution_clock::now();


    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration_generate = duration_cast<microseconds>(stop_generate - start_generate);

    cout << "Time taken for generating nodes, matrix, vector: "
         << duration_generate.count() << " microseconds" << endl;

    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration_solve = duration_cast<microseconds>(stop_solve - start_solve);

    cout << "Time taken for solving nodes, matrix, vector: "
         << duration_solve.count() << " microseconds" << endl;

//    check_solution(solution);
    write_solutions_to_file();

    return 0;
}

void check_solution(double *solution) {

    std::cout << "Comparing our results with sin(pi * x)" << std::endl;

    for(int i = 0; i < nodes_number; i++) {


        std::cout << solution[i] << " vs. " << a_sin_bx(1, PI, nodes[i + 1]) << std::endl;


        //checking vor an easier example
        // delete later

//        std::cout << "Comparing our results with 2,5*x^2" << std::endl;
//        std::cout << solution[i] << " vs. " << (0.5 * nodes[i + 1] - 0.5 * nodes[i + 1] * nodes[i + 1]) << std::endl;

    }

//    double solution_real[nodes_number];
//    for (int i = 1; i <= nodes_number; i++) {
//        for (int j = 1; j < nodes_number; j++) {
//        solution_real[i - 1] += solution[j - 1] * fi[i][j];
//            //std::cout << solution[j - 1] << " * " << fi[i][j] << std::endl;
//    }
//        //std::cout << solution_real[i - 1] <<std::endl;
//
//    }
}

void create_matrix() {
    create_nodes(nodes);

    //Randbedingungen setzten
    f_x_1d[0] = 0;
    f_x_1d[nodes_number + 1] = a_sin_bx(PI*PI , PI, (nodes[nodes_number + 1]+ nodes[nodes_number]) / 2);

    for (int i = 1; i <= nodes_number; i++) {



        //calculating the main function for vector b_i
        // f ((x_k-1 + x_k)/2)

        f_x_1d [i] = a_sin_bx(PI*PI , PI, (nodes[i]+ nodes[i - 1]) / 2);
        //f_x_1d [nodes_number] = 0;
        //delete later
//        f_x_1d [i] = 5;
//        f_x_1d [i] = 0.5 * (nodes[i + 1] - nodes[i + 1] * nodes[i + 1]);

        //calculating base function
        for (int j = 0; j <= nodes_number + 1; j ++) {
            fi[i][j]= calculate_fi((nodes[j] + nodes[j - 1]) / 2, nodes, i);
//            std::cout << "for base function " << i <<  " at node " << j << " for x = " << (nodes[j] + nodes[j - 1]) / 2 << " fi =  " << fi[i][j] << std::endl;
        }

        //calculating base function derivative
        //I switched i vs. j here to make it faster. This works because there are the same # functions and # points
        //for this loop i represents the number of the function, whereas j is the node number
        //in case those should not be in the same amount, this should be redesigned
        for (int j = 0; j <= nodes_number + 1; j ++) {
            dfi[i][j]= calculate_dfi((nodes[j] + nodes[j - 1]) / 2, nodes, i);
//            std::cout << "for base function " << i <<" at node " << j << " for x = " << (nodes[j] + nodes[j - 1]) / 2 << " dfi =  " << dfi[i][j] << std::endl;
        }

            //std::cout << " f at this point = " << f_x_1d[i]<< std::endl;

//        std::cout << "-----------------------------------------------" << std::endl;



    }

//    for (int i = 1; i <= nodes_number; i++) {
//        f_x_1d[i] = 1;
//        std::cout << "at k = " << i << " function value is " << f_x_1d[i] << std::endl;
//    }





//    for (int k = 1; k < nodes_number; k ++) {
//        std::cout << nodes[k] << std::endl;
//        for ( int z = 0; z < iter; z++)
////            std::cout << "-------------------------------------------------------"<< std::endl;
////
//            std::cout << f_x[k -1][z] << ", ";
//            std::cout << std::endl <<"-------------------------------------------------------"<< std::endl;
////            std::cout << "-------------------------------------------------------"<< std::endl;
//
//        std::cout << std::endl;
//    }


}




//function a * sin (b * x)
double a_sin_bx (double param1, double param2, double x) {

    double result = param1 * sin (param2 * x);
    //std::cout << param1 << " * " << "( sin (" << param2 << " * " << x << ") = " << result << std::endl;
    return result;

}



double calculate_dfi(double x, double *nodes, int i) {

//                                     __
//                                   (
//                                    \     (x - x_i-1)
//                                     \   -------------    for  x_i-1 <= x <= x_i     ==>     fi'(x) = 1 / (x_i - x_i-1)
//                                     |  (x_i - x_i-1)
//   for i = 1,..., N - 1, fi_i(x) =  /
//                                  <      (x_i+1 - x)
//                                   \     ----------      for   x_i <= x <= x_i+1     ==>     fi'(x) = -1 / (x_i+1 - x_i)
//                                   |    (x_i+1 - x_i)
//                                  /
//                                /       0                otherwise                   ==>     0
//                               (__
//
//        and                           __
//                                    (
//                                     \   (x - x_N-1)
//                                     |   -------------    for  x_N-1 <= x <= x_N     ==>     fi_N'(x) = -1 / (x_N - x_N-1)
//                                    /  (x_N - x_N-1)
//                       fi_N(x) =  <
//                                   \
//                                   |
//                                  /       0                otherwise                 ==>     0
//                                 (__

//
// x_i-1 <= x <= x_i?       fi_i = (x - x_i-1) / (x_i - x_i-1)
// x_i <= x <= x_i+1?       fi_i = (x_i+1 - x) / (x_i+1 - x_i)
// otherwise                fi_i = 0

    double dfi;
    if (nodes[i - 1] <= x && x <= nodes[i])
         dfi = 1 / (nodes[i] - nodes[i - 1]);
    else if (nodes[i] <= x && x <= nodes[i + 1] /* && i != nodes_number */)
         dfi = -1 / (nodes[i + 1] - nodes[i]);
    else dfi = 0;

//        std::cout << " node i - 1 = " << nodes[i - 1] << " node i = " << nodes[i] << " node i + 1 = " << nodes[i + 1] << std::endl;
//        std::cout << " at node " << i << " for x = " << x << " fi =  " << fi << std::endl;


        return dfi;


}


double calculate_fi(double x, double *nodes, int i) {

//                                     __
//                                   (
//                                    \     (x - x_i-1)
//                                     \   -------------    for  x_i-1 <= x <= x_i
//                                     |  (x_i - x_i-1)
//   for i = 1,..., N - 1, fi(x) =    /
//                                  <      (x_i+1 - x)
//                                   \     ----------      for   x_i <= x <= x_i+1
//                                   |    (x_i+1 - x_i)
//                                  /
//                                /       0                otherwise
//                               (__
//
//        and                           __
//                                    (
//                                     \   (x - x_N-1)
//                                     |   -------------    for  x_N-1 <= x <= x_N
//                                    /  (x_N - x_N-1)
//                         fi(x) =  <
//                                   \
//                                   |
//                                  /       0                otherwise
//                                 (__

//
// x_i-1 <= x <= x_i?       fi_i = (x - x_i-1) / (x_i - x_i-1)
// x_i <= x <= x_i+1?       fi_i = (x_i+1 - x) / (x_i+1 - x_i)
// otherwise                fi_i = 0

    double fi;
    if (nodes[i - 1] <= x && x <= nodes[i])
        fi = (x - nodes[i - 1]) / (nodes[i] - nodes[i - 1]);
    else if (nodes[i] <= x && x <= nodes[i + 1] /* && i != nodes_number */)
        fi = (nodes[i + 1] - x) / (nodes[i + 1] - nodes[i]);
    else fi = 0;

    //std::cout << " node i - 1 = " << nodes[i - 1] << " node i = " << nodes[i] << " node i + 1 = " << nodes[i + 1] << std::endl;
//        std::cout << " at node " << i << " for x = " << x << " fi =  " << fi << std::endl;


    return fi;


}

double fRand(double fMin, double fMax) {

    double f = (double) rand() / RAND_MAX;
    //std::cout << f_x;
    //std::cout << "fmin " <<fMin << " fmax " << fMax << " rand "<< f_x << " RAND_MAX " << RAND_MAX << std::endl;
    return fMin + f * (fMax - fMin);

}

void create_nodes (double *nodes) {
    nodes[0] = 0;
    nodes[nodes_number + 1] = 1;
//    std::cout << " x_0 = " << nodes[0] << std::endl;

    srand(time(NULL));
    rand();
    rand();
    rand();

    for (int i = 1; i <= nodes_number; i++) {

        //generating random nodes for base functions
        nodes[i] = fRand((double)fmax(nodes[i - 1], i - 1)/ (nodes_number + 1), (double) i/ (nodes_number + 1));
        //std::cout << "generated between " << fmax (nodes[i - 1], i - 1) << " and " << i + 1 << std::endl;
//        std::cout << " x_" << i << "= " << nodes[i] << std::endl;


        }
//    std::cout << " x_" << nodes_number+1 << "= " << nodes[nodes_number +1 ] << std::endl;

//    std::cout << " x_" << nodes_number << "= " << nodes[nodes_number] << std::endl;
//    for (int i = 1; i <= nodes_number; i++ ) {
//
//        nodes[i] = (double)i / (nodes_number + 1);
//        std::cout << " x_" << i << "= " << nodes[i] << std::endl;
//
//    }
//    std::cout << " x_" << nodes_number + 1 << "= " << nodes[nodes_number + 1] << std::endl;


}

void create_matrix_b(double *b, double *f_1d) {

    //---------------------------------------------------------------------------------------------------------
    //             N
    //            ___
    //      b_i = \                            x_k-1 + x_k                 x_k-1 + x_k
    //            /     (x_k - x_k-1)  *  f ( ------------- )  *   fi_i ( ------------- )
    //           ----                             2                            2
    //           k=1
    //
    //      for i = 1,..., node_number.
    //
    //---------------------------------------------------------------------------------------------------------



    for (int i = 1; i <= nodes_number; i++) {
        for (int k = 0; k <= nodes_number; k++) {
            b[i - 1] += (nodes[k + 1] - nodes[k]) * f_x_1d[k + 1] * fi[i][k + 1];
//            std::cout << "(" << nodes[k + 1] << " - " << nodes[k] << ") * " << f_x_1d[k + 1]  << " * " << fi[i][k + 1] << std::endl;
        }
//        std::cout << " b[" << i - 1 << "] = " << b[i - 1] << std::endl;
    }

}

void create_matrix_K() {

    //---------------------------------------------------------------------------------------------------------
    //
    //             _N__
    //      K_ij = \                               x_k-1 + x_k                  x_k-1 + x_k
    //            /     (x_k - x_k-1)  *  fi_i' ( ------------- )  *   fi_j' ( ------------- )
    //           ----                                  2                            2
    //           k=1
    //
    //      for i = 1,..., node_number.
    //
    //---------------------------------------------------------------------------------------------------------

//    std::cout << "-----------------------------------------------" << std::endl;

    for (int i = 1; i <= nodes_number; i++) {
        for (int j = 1; j <= nodes_number; j++) {
            for (int k = 0; k <= nodes_number; k++) {
                K[i - 1][j - 1] += (nodes[k + 1] - nodes[k]) * dfi[i][k + 1] * dfi[j][k + 1];
//                std::cout << "(" << "nodes[" << k << "]" << " - " << "nodes["<<k-1<< "]" << ") * " << "dfi[" << i << "][" <<  k << "]"  << " *  dfi[" << j << "][" << k << "]" << std::endl;
//                std::cout << "(" << nodes[k] << " - " << nodes[k-1] << ") * " << dfi[i][k]  << " * " << dfi[j][k] << std::endl;
            }
//            std::cout << " ------------" << i << "  " << j << "------------------------- " << std::endl;
//            std::cout << " K[" << i << "][" << j << "] = " << K[i - 1 ][j - 1];
//            std::cout <<"\t"<< K[i - 1][j - 1];

        }
//        std::cout << std::endl;
    }

//    std::cout << "-----------------------------------------------" << std::endl;

}


double * solve_matrix () {

//        timestamp ( );
//        cout << "\n";
//        cout << "SOLVE_TEST\n";
//        cout << "  C++ version\n";
//        cout << "  Test the SOLVE library.\n";



    //****************************************************************************80
//
//  Purpose:
//
//    TEST01 demonstrates how a 3X3 linear system can be set up and solved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 May 2014
//
//  Author:
//
//    John Burkardt
//

        double **a;
        int n;
        double *x;
//
//  Define the array size.
//
        n = nodes_number;
//
//  Create the array that will contain the matrix.
//
        a = r8rmat_new ( n, n );
//
//  Set the array values.
//

    for(int i = 0; i < nodes_number; i++) {
        a[i] = K[i];
//        for (int j = 0; j < nodes_number; j++)
//            std:: cout << a[i][j] << " = " << K[i][j] <<std::endl;
    }

//    for(int i = 1; i <= nodes_number; i++) {
//        b_solver[i - 1] = b[i];
//        for (int j = 0; j < nodes_number; j++)
//            std:: cout << a[i][j] << " = " << K[i][j] <<std::endl;
//    }

//
//  Create the right hand side.
//
//        b = new double[n];
////
////  Set the right hand side values.
////
//        b[0] = 14;
//        b[1] = 32;
//        b[2] = 23;
//
//  Request the solution of A*x=b.

//        for (int j = 0; j < nodes_number; j++)
//            std::cout << b[j] << std::endl;
//            std::cout << " b[" << j << "] = " << b[j] << std::endl;
//
        x = r8rmat_fs_new ( n, a, b );

//        r8vec_print ( n, x, "  Solution:" );
//
//  Free memory.
//
//        r8rmat_delete ( n, n, a );
//        delete [] b;
//        delete [] x;


        //
//  Terminate.
//
//        cout << "\n";
//        cout << "SOLVE_TEST\n";
//        cout << "  Normal end of execution.\n";
//        cout << "\n";
//        timestamp ( );



        return x;
}

void write_solutions_to_file() {

    //
//  Open data file, and write solutions as they are computed.
//
    data_unit.open ( data_filename.c_str ( ) );


//    data_unit << "x " << "y " << "velocity" <<"\n";
//    data_unit << "#coordflag xv" <<"\n";
    for (int j = 0; j <= nodes_number; j++ )
    {
        data_unit << nodes[j + 1]
                  << " " << solution[j]
                  << "\n";
    }
    data_unit << "\n";

}