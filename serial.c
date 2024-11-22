#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_N 10000  // Increased maximum value for N

void lu_decomposition_doolittle(int A[MAX_N][MAX_N], float L[MAX_N][MAX_N], float U[MAX_N][MAX_N], int N) {
    int i, j, k;

    for (i = 0; i < N; i++) {
        L[i][i] = 1.0;

        // Compute U[i][j]
        for (j = i; j < N; j++) {
            float sum = 0.0;
            for (k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = A[i][j] - sum;
        }

        // Compute L[j][i]
        for (j = i + 1; j < N; j++) {
            float sum = 0.0;
            for (k = 0; k < i; k++) {
                sum += L[j][k] * U[k][i];
            }
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }
}

int main() {
    int N, num_threads, i, j;
    int A[MAX_N][MAX_N];
    float L[MAX_N][MAX_N] = {0}, U[MAX_N][MAX_N] = {0};

    // Input matrix size N
    printf("Enter the size of the matrix (N): ");
    scanf("%d", &N);

    if (N <= 0 || N > MAX_N) {
        printf("Invalid matrix size. Exiting...\n");
        return 1;
    }

    // Seed for random number generation
    srand(time(NULL));

    // Generate random matrix A between -100 and 100
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[i][j] = rand() % 201 - 100;
        }
    }

    // Display Matrix A
    printf("\nMatrix A:\n");
  /*  for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }*/

    // Measure time taken for LU decomposition
    double start_time = (double)clock() / CLOCKS_PER_SEC;
    lu_decomposition_doolittle(A, L, U, N);
    double end_time = (double)clock() / CLOCKS_PER_SEC;

  /*  // Displaying L, U matrices
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
    }*/

    // Display time taken
    printf("\nTime taken: %.6f seconds\n", end_time - start_time);

    return 0;
}

