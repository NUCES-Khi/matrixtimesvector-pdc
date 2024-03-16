# PDC Assignments and Projects
## Team Members
|std_id|Name|
|--------|-|
|k21-3054|Khadija Ibrahim|
|k21-4735|Nujood Idrees|
|k21-4702|Abdul Basit|

## Assingment 1 ##
Status: **completed**
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

// Function to perform sequential matrix-vector multiplication
void matrixVectorMultiplicationSequential(double **matrix, double *vector, double *result, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

// Function to perform matrix-vector multiplication using OpenMP
void matrixVectorMultiplicationOpenMP(double **matrix, double *vector, double *result, int rows, int cols) {
    #pragma omp parallel for
    for (int i = 0; i < rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

// Function to perform matrix-vector multiplication using MPI
void matrixVectorMultiplicationMPI(double **matrix, double *vector, double *result, int rows, int cols, int rank, int num_procs) {
    int chunk_size = rows / num_procs;
    int start_row = rank * chunk_size;
    int end_row = (rank == num_procs - 1) ? rows : (rank + 1) * chunk_size;

    for (int i = start_row; i < end_row; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <matrix_size> <vector_size>\n", argv[0]);
        return 1;
    }

    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);
    int rank, num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Allocate memory for matrix, vector, and result
    double **matrix = (double **)malloc(rows * sizeof(double *));
    double *vector = (double *)malloc(cols * sizeof(double));
    double *result = (double *)malloc(rows * sizeof(double));

    // Fill matrix and vector with random values (only rank 0)
    if (rank == 0) {
        for (int i = 0; i < rows; i++) {
            matrix[i] = (double *)malloc(cols * sizeof(double));
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = (double)rand() / RAND_MAX;
            }
        }

        for (int i = 0; i < cols; i++) {
            vector[i] = (double)rand() / RAND_MAX;
        }
    }

    // Broadcast vector to all processes
    MPI_Bcast(vector, cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Perform matrix-vector multiplication
    double start_time, end_time, execution_time;
    start_time = MPI_Wtime();

    // Choose implementation based on rank
    if (num_procs == 1 || cols <= 1000) {
        // Use sequential implementation
        matrixVectorMultiplicationSequential(matrix, vector, result, rows, cols);
    } else if (cols <= 10000) {
        // Use OpenMP implementation
        matrixVectorMultiplicationOpenMP(matrix, vector, result, rows, cols);
    } else {
        // Use MPI implementation
        matrixVectorMultiplicationMPI(matrix, vector, result, rows, cols, rank, num_procs);
    }

    end_time = MPI_Wtime();
    execution_time = end_time - start_time;

    // Gather results from all processes to rank 0
    MPI_Gather(result, rows / num_procs, MPI_DOUBLE, result, rows / num_procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Print result and execution time (only rank 0)
    if (rank == 0) {
        printf("Result vector:\n");
        for (int i = 0; i < rows; i++) {
            printf("%.2f ", result[i]);
        }
        printf("\n");

        printf("Execution time: %.6f seconds\n", execution_time);
    }

    // Deallocate memory
    if (rank == 0) {
        for (int i = 0; i < rows; i++) {
            free(matrix[i]);
        }
        free(matrix);
        free(vector);
    }
    free(result);

    MPI_Finalize();

    return 0;
}

**Batch Script**
@echo off

REM Define input sizes
set input_sizes=64 128 256 512 1024 2048 4096 8192 16384

REM Number of runs for each input size
set num_runs=10

REM Output file name
set output_file=benchmark_results.csv

REM Clear existing output file
echo test S.no, file, input size, time taken, average so far > %output_file%

REM Loop through each input size
for %%s in (%input_sizes%) do (
    REM Loop through each program
    for %%p in ("sequential" "openmp" "mpi") do (
        REM Run program num_runs times and record execution times
        for /l %%i in (1,1,%num_runs%) do (
            echo Running %%p program for input size %%s, Run %%i
            (time (%%p.exe %%s %%s 2>&1)) >> temp.txt
            REM Extract execution time from temp.txt
            for /f "tokens=2" %%t in ('findstr "real" temp.txt') do (
                REM Write execution time to output file
                echo %%p.exe, %%s, %%t >> %output_file%
            )
        )
        REM Delete temporary file
        del temp.txt
    )
)

echo Benchmarking complete!

