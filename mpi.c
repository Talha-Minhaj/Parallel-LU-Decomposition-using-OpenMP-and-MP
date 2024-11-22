#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define MAX_N 10000  // Maximum value for N

void lu_decomposition_doolittle(int A[MAX_N][MAX_N], float L[MAX_N][MAX_N], float U[MAX_N][MAX_N], int N, int rank, int num_processes) {
    int i, j, k;

    for (i = 0; i < N; i++) {
        L[i][i] = 1;

        // Distribute the row work across processes
        for (j = i; j < N; j++) {
            U[i][j] = A[i][j];
            for (k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }

        // Distribute the column work across processes
        for (j = i + 1; j < N; j++) {
            L[j][i] = A[j][i];
            for (k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
}

int main(int argc, char *argv[]) {
    int N, i, j;
    int A[MAX_N][MAX_N];
    int rank, num_processes;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    // Seed for random number generation
    srand(time(NULL) + rank);  // Ensure different seed for each process

    // Input matrix size N
    if (rank == 0) {
        printf("Enter the size of the matrix (N): ");
        scanf("%d", &N);

        if (N <= 0 || N > MAX_N) {
            printf("Invalid matrix size. Exiting...\n");
            MPI_Finalize();
            return 1;
        }
    }

    // Broadcast matrix size to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Generate random matrix A between -100 and 100 on process 0
    if (rank == 0) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                A[i][j] = rand() % 201 - 100;
            }
        }
    }

    // Broadcast the matrix A to all processes
    MPI_Bcast(A, N * N, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate L and U matrices
    float L[MAX_N][MAX_N] = {0}, U[MAX_N][MAX_N] = {0};

    // Measure time taken for LU decomposition
    double start_time = MPI_Wtime();
    lu_decomposition_doolittle(A, L, U, N, rank, num_processes);
    double end_time = MPI_Wtime();

    // Display Matrix A on rank 0
    if (rank == 0) {
   /*     printf("\nMatrix A:\n");
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%d ", A[i][j]);
            }
            printf("\n");
        }

        // Displaying L and U matrices on rank 0
        printf("\nMatrix L:\n");
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%.2f ", L[i][j]);
            }
            printf("\n");
        }

        printf("\nMatrix U:\n");
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%.2f ", U[i][j]);
            }
            printf("\n");
        }
*/
        // Display time taken
        printf("\nTime taken: %.6f seconds\n", end_time - start_time);
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}

