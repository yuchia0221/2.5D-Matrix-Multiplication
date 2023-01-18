#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

int c;
const int X_AXIS = 0;
const int Y_AXIS = 1;
const int FRONT_PLANE = 0;

double l2_norm(double);
bool check_inputs(int, int, int);

void print_matrix(double *&, int);
void initialize_matrices(double *&, double *&, double *&);
void block_multiplication(double *&, double *&, double *&, int);
void matrix_multiplication(double *&, double *&, double *&, int);
void verify_multiplication_result(double *&, double *&, double *&, int);
void initialize_buffers(double *&, double *&, double *&, double *&, int);
void restore_matrix_order(double *&, double *&, double *&, double *&, double *&, double *&, int);

void print_matrix(double *&matrix, int size)
{
    // Debug only: print the whole matrix
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
            printf("%.2lf ", matrix[i * size + j]);
        cout << endl;
    }
    cout << endl;
}

double l2_norm(double value)
{
    // Return the value of l2 norm
    return sqrt(pow(value, 2));
}

bool check_inputs(int processor_num, int c, int matrix_size)
{
    double rounds = sqrt(processor_num / pow(c, 3));
    double block_size = matrix_size / sqrt(processor_num / c);
    return ceilf(rounds) == rounds && ceilf(block_size) == block_size;
}

void initialize_matrices(double *&A, double *&B, double *&C, int size)
{
    // Initialize matrices (A, B)
    A = new double[size * size];
    B = new double[size * size];
    C = new double[size * size];

    // Set random seed
    srand48(221);

    // Randomly generate double percision matrix with 0 <= value <= 10
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            A[i * size + j] = (int)(drand48() * 10);
            B[i * size + j] = (int)(drand48() * 10);
            C[i * size + j] = 0;
        }
    }
}

void initialize_buffers(double *&buffer1, double *&buffer2, double *&buffer3, double *&buffer4, int size)
{
    // Initialize Buffer Array for MPI commuunication
    buffer1 = new double[size]{};
    buffer2 = new double[size]{};
    buffer3 = new double[size]{};
    buffer4 = new double[size]{};
}

void block_multiplication(double *&A, double *&B, double *&C, int block_size)
{
    // Cache-friendly matrix multiplication of two blocks
    for (int k = 0; k < block_size; k++)
    {
        for (int i = 0; i < block_size; i++)
        {
            double r = A[i * block_size + k];
            for (int j = 0; j < block_size; j++)
            {
                C[i * block_size + j] += r * B[k * block_size + j];
            }
        }
    }
}

void matrix_multiplication(double *&A, double *&B, double *&C, int matrix_size)
{
    double start, end;
    int ndims = 3, reorder = 1;
    int world_rank, world_size, coords[3];

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int block_size = matrix_size / sqrt(world_size / c);
    int element_number = block_size * block_size;
    int grid_dim = matrix_size / block_size;
    int dims[3] = {grid_dim, grid_dim, c}, periods[3] = {1, 1, 0};

    // Initailize buffer array for MPI communication
    double *A_buffer, *B_buffer, *C_buffer, *C_result;
    initialize_buffers(A_buffer, B_buffer, C_buffer, C_result, element_number);

    // Set up MPI communicators
    MPI_Comm comm_cart;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);
    MPI_Cart_coords(comm_cart, world_rank, ndims, coords);

    start = MPI_Wtime();

    // Scatter A and B matrix to the front plane
    MPI_Comm comm_xy_cart;
    MPI_Cart_sub(comm_cart, new int[3]{1, 1, 0}, &comm_xy_cart);
    if (coords[2] == FRONT_PLANE)
    {
        MPI_Scatter(A, element_number, MPI_DOUBLE, A_buffer, element_number, MPI_DOUBLE, 0, comm_xy_cart);
        MPI_Scatter(B, element_number, MPI_DOUBLE, B_buffer, element_number, MPI_DOUBLE, 0, comm_xy_cart);
    }

    // Broadcast A and B among the Z-axis
    MPI_Comm comm_z_cart;
    MPI_Cart_sub(comm_cart, new int[3]{0, 0, 1}, &comm_z_cart);
    MPI_Bcast(A_buffer, element_number, MPI_DOUBLE, 0, comm_z_cart);
    MPI_Bcast(B_buffer, element_number, MPI_DOUBLE, 0, comm_z_cart);

    // Initial Alignment for performing Cannon Algorithm
    int rounds = grid_dim / c;
    int r = (coords[1] + coords[0] - coords[2] * rounds) % grid_dim;
    int s = (coords[1] - coords[0] + coords[2] * rounds) % grid_dim;
    int s1 = (coords[0] - coords[1] + coords[2] * rounds) % grid_dim;
    int A_recieve, A_send, B_recieve, B_send;
    MPI_Cart_rank(comm_cart, new int[3]{coords[0], s, coords[2]}, &A_send);
    MPI_Cart_rank(comm_cart, new int[3]{coords[0], r, coords[2]}, &A_recieve);
    MPI_Cart_rank(comm_cart, new int[3]{s1, coords[1], coords[2]}, &B_send);
    MPI_Cart_rank(comm_cart, new int[3]{r, coords[1], coords[2]}, &B_recieve);

    // Send and Receive blocks of A and B matrix from other processes
    MPI_Status A_status, B_status;
    MPI_Sendrecv_replace(A_buffer, element_number, MPI_DOUBLE, A_send, 1, A_recieve, 1, comm_cart, &A_status);
    MPI_Sendrecv_replace(B_buffer, element_number, MPI_DOUBLE, B_send, 1, B_recieve, 1, comm_cart, &B_status);

    // Multiply two alignmnet blocks
    block_multiplication(A_buffer, B_buffer, C_buffer, block_size);

    // Perform Cannon's Algorithm
    MPI_Cart_shift(comm_cart, Y_AXIS, 1, &A_recieve, &A_send);
    MPI_Cart_shift(comm_cart, X_AXIS, 1, &B_recieve, &B_send);
    for (int i = 0; i < rounds - 1; i++)
    {
        MPI_Sendrecv_replace(A_buffer, element_number, MPI_DOUBLE, A_send, 1, A_recieve, 1, comm_cart, &A_status);
        MPI_Sendrecv_replace(B_buffer, element_number, MPI_DOUBLE, B_send, 1, B_recieve, 1, comm_cart, &B_status);
        block_multiplication(A_buffer, B_buffer, C_buffer, block_size);
    }

    // Gather matrix multiplication result
    MPI_Reduce(C_buffer, C_result, element_number, MPI_DOUBLE, MPI_SUM, 0, comm_z_cart);
    if (coords[2] == FRONT_PLANE)
        MPI_Gather(C_result, element_number, MPI_DOUBLE, C, element_number, MPI_DOUBLE, 0, comm_xy_cart);

    // Print Execution Time
    if (world_rank == 0)
    {
        end = MPI_Wtime();
        printf("2.5D Matrix Multiplication finished in: %.4f seconds\n", end - start);
    }

    // Free Memory
    delete[] A_buffer;
    delete[] B_buffer;
    delete[] C_buffer;
    delete[] C_result;
};

void restore_matrix_order(double *&A, double *&B, double *&C, double *&A_restore, double *&B_restore, double *&C_restore, int matrix_size)
{
    // Restore matrix to its orignal order
    A_restore = new double[matrix_size * matrix_size]{};
    B_restore = new double[matrix_size * matrix_size]{};
    C_restore = new double[matrix_size * matrix_size]{};

    int processors;
    MPI_Comm_size(MPI_COMM_WORLD, &processors);

    int counter = 0, rounds = 0, row_index = 0;
    int block_dim = matrix_size / sqrt(processors / c);
    int grid_dim = matrix_size / block_dim;

    for (int i = 0; i < processors / c; i++)
    {
        rounds = i % grid_dim;
        row_index = (i % (processors / c)) / grid_dim * block_dim;
        for (int k = 0; k < block_dim; k++)
        {
            for (int j = 0; j < block_dim; j++)
            {
                int index = row_index * matrix_size + j + rounds * block_dim;
                A_restore[index] = A[counter];
                B_restore[index] = B[counter];
                C_restore[index] = C[counter];
                ++counter;
            }
            ++row_index;
        }
    }
}

void verify_multiplication_result(double *&A, double *&B, double *&C, int matrix_size)
{
    // Verify for the result of 2.5D Matrix Multiplication
    double *A_restore, *B_restore, *C_restore;
    restore_matrix_order(A, B, C, A_restore, B_restore, C_restore, matrix_size);

    double residual;
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            for (int k = 0; k < matrix_size; k++)
            {
                C_restore[i * matrix_size + j] -= A_restore[i * matrix_size + k] * B_restore[k * matrix_size + j];
            }
            residual += l2_norm(C_restore[i * matrix_size + j]);
        }
    }
    printf("The sum of Euclidean norms is: %.2e\n", residual);

    // Free Memory
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] A_restore;
    delete[] B_restore;
    delete[] C_restore;
}

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);
    double *A, *B, *C;
    int world_rank, world_size, matrix_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    matrix_size = atoi(argv[1]), c = atoi(argv[2]); // Get input argument

    if (!check_inputs(world_size, c, matrix_size))
    {
        if (world_rank == 0)
            cout << "Input is not valid: sqrt(p/c) has to be a multiple of c" << endl;
        exit(1);
    }

    if (world_rank == 0)
    {
        printf("Using %d processes to run 2.5D matrix multiplication with arrangement %d*%d*%d (c*x*y)\n",
               world_size, c, (int)sqrt(world_size / c), (int)sqrt(world_size / c));
        initialize_matrices(A, B, C, matrix_size);
    }

    matrix_multiplication(A, B, C, matrix_size);

    if (world_rank == 0)
        verify_multiplication_result(A, B, C, matrix_size);

    MPI_Finalize();

    return 0;
}