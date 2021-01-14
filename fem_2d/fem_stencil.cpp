#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

// solving Uxx+Uyy = F
// [K] * [solution] = [b]
// we are considering Dirichlet BC for simplicity
// We use sparse striffness matrix.

using namespace std;
using namespace chrono;

#define PI 3.14159265358979323846
constexpr int nodes_number_x{8};
constexpr int nodes_number_y{8};
constexpr int nodes_number{nodes_number_x * nodes_number_y};

// #define DBG

// BC: 
/* mesh
T - T - T - T
|   |   |   |
L - x - x - R
|   |   |   |  
L - x - x - R
|   |   |   |
L - B - B - R 
*/

/* memory layout for the sparse stiffness matrix:
size of sparse matrix = 9 * nodes_number
9-point stencil -> we have 9 elements per row
1  1  1  1  1  1  1  1  1 | 2  2  2  2  2  2  2  2  2 | 3  3  3  3  3  3  3  3  3  | ...
-1 -1 -1 -1 8 -1 -1 -1 -1 | -1 -1 -1 -1 8 -1 -1 -1 -1 | -1 -1 -1 -1 8 -1 -1 -1 -1 | ... 
* */

//https://www5.in.tum.de/lehre/vorlesungen/sci_comp/ws15/pde_fem-2x4.pdf

/* nodes
x - x - x - x
|   |   |   |
x - x - x - x
|   |   |   |  
x - x - x - x
|   |   |   |
x - x - x - x 
*/

/* F values 
x - x - x - x
|   |   |   |
x - x - x - x
|   |   |   |  
x - x - x - x
|   |   |   |
x - x - x - x 
*/

/* indirection_m: offsets used in the sparse representation of matrix
*/

/* indirection_v: indirections used in the sparse representation to get columns
*/

// helper functions
void create_grid();
double a_sin_bxy(double param1, double param2, double x, double y);
void create_fxy(double *f);

// FEM function
void create_stencil(double scalor);
void create_indirections();
void create_matrix_b(double *b, double *f, double scalor);

// write the results
void write_solutions_to_file();

ofstream command_unit;
string data_filename = "fem_data.csv";
ofstream data_unit;

// grid buffers
double nodes_x[nodes_number];
double nodes_y[nodes_number];

// buffer for the source function
double f_xy[nodes_number];

// buffer for b vector
double b[nodes_number];

// stiffness matrix
//double stiffness_matrix[nodes_number * 9]{};
// stencil tensor
double stencil[9]{};

// indirections
int indirections[nodes_number * 9]{};

double solution_1[nodes_number]{};
double solution_2[nodes_number]{};

int main()
{
    // Get starting timepoint
    auto start_generate = high_resolution_clock::now();
    
    create_grid();

    create_fxy(f_xy);

    int scalor_b = 1.0;
    create_matrix_b(b, f_xy, scalor_b);

    int scalor_K = 1.0;
    // create_sparse_stiffness_matrix(scalor_K);
    create_stencil(scalor_K);
    create_indirections();

    // cout << "Comparing our results with sin(pi * x)" << endl;

#ifdef DBG
    cout << "nodex_x: " << endl;
    for (int i = 0; i < nodes_number; i++)
    {
        cout << nodes_x[i] << " " << endl;
    }
    cout << endl;

    cout << "nodex_y: " << endl;
    for (int i = 0; i < nodes_number; i++)
    {
        cout << nodes_y[i] << " " << endl;
    }
    cout << endl;

    cout << "stencil: " << endl;
    for (int i = 0; i < 9; i++)
    {
        cout << stencil[i] << " " << endl;
    }
    cout << endl;

    cout << "indirections: " << endl;
    for (int i = 0; i < nodes_number * 9; i++)
    {
        if (i % 9 == 0)
        {
            cout << "row: ith: " << i / 9 << endl;
        }
        cout << indirections[i] << " " << endl;
    }
    cout << endl;

    cout << "f: " << endl;
    for (int i = 0; i < nodes_number; i++)
    {
        cout << f_xy[i] << " " << endl;
    }
    cout << endl;

    cout << "b: " << endl;
    for (int i = 0; i < nodes_number; i++)
    {
        cout << b[i] << " " << endl;
    }
    cout << endl;
#endif

    // Get ending timepoint
    auto stop_generate = high_resolution_clock::now();

    // Get starting timepoint
    auto start_solve = high_resolution_clock::now();
    // solution = solve_matrix();
    for (int i = 0; i < nodes_number; i++)
    {
        solution_1[i] = 0.0;
    }

    int num_iters = 100;

    // iterations of jacobi
    double *s1 = solution_1;
    double *s2 = solution_2;
    double *tmp;
    for (size_t it = 0; it < num_iters; it++)
    {
        // iterate through all the rows of the matrix
        /*
        for (size_t i = 0; i < nodes_number; i++)
        {
            s2[i] = 1.0 / stencil[4] * (b[i] - 
                (stencil[0] * s1[indirections[9 * i + 0]] + 
                 stencil[1] * s1[indirections[9 * i + 1]] + 
                 stencil[2] * s1[indirections[9 * i + 2]] + 
                 stencil[3] * s1[indirections[9 * i + 3]] + 
                 stencil[5] * s1[indirections[9 * i + 5]] + 
                 stencil[6] * s1[indirections[9 * i + 6]] + 
                 stencil[7] * s1[indirections[9 * i + 7]] + 
                 stencil[8] * s1[indirections[9 * i + 8]] ) );
        }
        */

        // simple jacobi solver
        // note Dirichlet BC is implied
        for (size_t idx_i = 1; idx_i < nodes_number_x - 1; idx_i++)
        {
            for (size_t idx_j = 1; idx_j < nodes_number_y - 1; idx_j++)
            {
                int i = idx_i * nodes_number_y + idx_j;
                s2[i] = 1.0 / stencil[4] * (b[i] - (stencil[0] * s1[indirections[9 * i + 0]] + stencil[1] * s1[indirections[9 * i + 1]] + stencil[2] * s1[indirections[9 * i + 2]] + stencil[3] * s1[indirections[9 * i + 3]] + stencil[5] * s1[indirections[9 * i + 5]] + stencil[6] * s1[indirections[9 * i + 6]] + stencil[7] * s1[indirections[9 * i + 7]] + stencil[8] * s1[indirections[9 * i + 8]]));
            }
        }

        // swap the buffers
        tmp = s1;
        s1 = s2;
        s2 = tmp;
    }

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

#ifdef DBG
    cout << "solution: " << endl;
    for (int i = 0; i < nodes_number; i++)
    {
        cout << solution_1[i] << " " << endl;
    }
    cout << endl;
#endif

    write_solutions_to_file();

    return EXIT_SUCCESS;
}

void create_fxy(double *f)
{
    for (int i = 0; i < nodes_number; i++)
    {
        f_xy[i] = a_sin_bxy(1.0, PI, nodes_x[i], nodes_y[i]);
    }
}

// function a * sin (b * x) * sin (b * y)
double a_sin_bxy(double param1, double param2, double x, double y)
{
    double result = param1 * sin(param2 * x) * sin(param2 * y);
#ifdef DBG
    cout << param1 << " * "
         << "( sin (" << param2 << " * " << x << ") = " << result << endl;
#endif
    return result;
}

void create_grid()
{
    for (int i = 0; i < nodes_number_x; i++)
    {
        for (int j = 0; j < nodes_number_x; j++)
        {
            nodes_x[i * nodes_number_y + j] = (double)j / (nodes_number_x - 1);
            nodes_y[i * nodes_number_y + j] = (double)i / (nodes_number_y - 1);
#ifdef DBG
            cout << "nodes_x[" << i << "*" << nodes_number_y << "+ " << j << "] = " << j << endl; //(double)i / (nodes_number_x - 1);
            cout << "nodes_y[" << i << "*" << nodes_number_y << "+ " << j << "] = " << i << endl; //(double)i / (nodes_number_x - 1);
#endif
        }
    }
}

void create_matrix_b(double *b, double *f, double scalor)
{
    for (int i = 0; i < nodes_number; i++)
    {
        b[i] = scalor * f[i];
    }
}

void create_stencil(double scalor)
{
    for (int i = 0; i < nodes_number; i++)
    {
        stencil[0] = -1 * scalor;
        stencil[1] = -1 * scalor;
        stencil[2] = -1 * scalor;
        stencil[3] = -1 * scalor;
        stencil[4] = 8 * scalor;
        stencil[5] = -1 * scalor;
        stencil[6] = -1 * scalor;
        stencil[7] = -1 * scalor;
        stencil[8] = -1 * scalor;
    }
}

void create_indirections()
{
    for (int i = 0; i < nodes_number; i++)
    {
        indirections[i * 9 + 0] = i + -1 - nodes_number_x;
        indirections[i * 9 + 1] = i + 0 - nodes_number_x;
        indirections[i * 9 + 2] = i + 1 - nodes_number_x;
        indirections[i * 9 + 3] = i + -1;
        indirections[i * 9 + 4] = i + 0;
        indirections[i * 9 + 5] = i + 1;
        indirections[i * 9 + 6] = i + -1 + nodes_number_x;
        indirections[i * 9 + 7] = i + 0 + nodes_number_x;
        indirections[i * 9 + 8] = i + 1 + nodes_number_x;
    }
}

void write_solutions_to_file()
{
    data_unit.open(data_filename.c_str());

    for (int i = 0; i < nodes_number; i++)
    {
        data_unit << nodes_x[i] << " " << nodes_y[i] << " " << solution_1[i] << "\n";
    }

    data_unit.flush();
    data_unit.close();
}