#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void print_matrix_to_file(double *matrix, int n, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%f ", matrix[i * n + j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void seq_mat_mult(double* mat1, double* mat2, double* mat3, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                mat3[i * n + j] += mat1[i * n + k] * mat2[k * n + j];
}

void transpose_matrix(int N, double *original_matrix, double *transposed_matrix)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            transposed_matrix[j * N + i] = original_matrix[i * N + j];
        }
    }
}

double multiply_row_and_col(int N, double *row, double *col)
{
    double total = 0;
    for (int i = 0; i < N; i++)
    {
        total += row[i] * col[i];
    }
    return total;
}

void check_correctness(double *result_seq, double *result_parallel, int n)
{
/*    FILE *file_seq, *file_parallel;
    char filename_seq[50], filename_parallel[50];

    // Print result_seq to file
    sprintf(filename_seq, "result_seq.txt");
    file_seq = fopen(filename_seq, "w");
    for (int i = 0; i < n; i++)
    {
         for (int j = 0; j < n; j++)
         {
              fprintf(file_seq, "%f ", result_seq[i * n + j]);
         }
         fprintf(file_seq, "\n");
    }
    fclose(file_seq);

    // Print result_parallel to file
    sprintf(filename_parallel, "result_parallel.txt");
    file_parallel = fopen(filename_parallel, "w");
    for (int i = 0; i < n; i++)
    {
         for (int j = 0; j < n; j++)
         {
              fprintf(file_parallel, "%f ", result_parallel[i * n + j]);
         }
         fprintf(file_parallel, "\n");
    }
    fclose(file_parallel);
*/
    int check_positive = 0;
    int check_count = 0;

    // Loop to check 10 values
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            // Compare the absolute difference of corresponding elements
             double diff = fabs(result_seq[i * n + j] - result_parallel[i * n + j]);

            // Print information about the comparison
            printf("Values at (%d, %d): Result_seq = %e, Result_parallel = %e, Diff = %e\n",i, j, result_seq[i * n + j], result_parallel[i * n + j], diff);

            // Check if the difference is within the threshold
            if (diff < 1e-10) {
                printf("Values at (%d, %d) are close to 0. Diff: %e\n", i, j, diff);
                check_positive++;
            }
            else{
                printf("Values at (%d, %d) are NOT close to 0. Diff: %e\n", i, j, diff);
            }
	check_count++;    
        }
    }

    // Print the total number of values checked
    printf("Total values checked: %d, close to 0 (positive check) are: %d\n", check_count, check_positive);
}

void parallel_mat_mult(double *mat1, double *mat2, double *mat3, int n, int size)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Calculate local matrix size
    int local_n = n / size;

    // Allocate memory for local matrices
    double *local_A = (double *)malloc(local_n * n * sizeof(double));
    double *local_B = (double *)malloc(local_n * n * sizeof(double));
    double *local_C = (double *)malloc(local_n * n * sizeof(double));


//    printf("Number check!!!!!!!!!! n: %d, local_n: %d, size: %d, rank: %d \n", n, local_n, size, rank);
    // Distribute data using MPI_Scatter
    MPI_Scatter(mat1, local_n * n, MPI_DOUBLE, local_A, local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(mat2, local_n * n, MPI_DOUBLE, local_B, local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Perform the Ring algorithm for matrix multiplication
    for (int step = 0; step < size; step++)
    {
	int dest = (rank - 1 + size) % size;
        int src = (rank + 1) % size;
	double sum = 0.0;

        // Perform matrix multiplication for this step
        for (int i = 0; i < local_n; i++)
        {    
	     for(int j = 0; j < local_n; j++)
	     {	     
                 double sum = 0.0;
                 for (int k = 0; k < n; k++)
                 {
                     sum += local_A[i * n + k] * local_B[j * n + k];
                 }
                 //local_C[(((local_n * (rank + step)) + i) % n) + i * n + (((local_n * (rank + step)) + j) % n)] = sum;
		//local_C[((local_n * (rank + step)) + i) * n + i + ((local_n * (rank + step)) + j)] = sum;
		 local_C[(i * n + j) + (rank * local_n + step * local_n) % n] = sum;
	     }
        }


        // Rotate the rows of local_B among processes
        if (step < size - 1)
        {
            MPI_Sendrecv_replace(local_B, local_n * n, MPI_DOUBLE, dest, 0, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
	MPI_Barrier(MPI_COMM_WORLD);
    }
   //the following is for debuggung
   /* 
    // Write local_C to a file unique to each processor
    char filename[50];
    sprintf(filename, "local_C_proc_%d.txt", rank);
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n", filename);
        return;
    }

    fprintf(file, "Processor %d local_C matrix:\n", rank);
    for (int i = 0; i < local_n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%f ", local_C[i * local_n + j]);
        }
        fprintf(file, "\n");
    }
    fclose(file); */

    // Gather results using MPI_Gather
    MPI_Gather(local_C, local_n * n, MPI_DOUBLE, mat3, local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Clean up
    free(local_A);
    free(local_B);
    free(local_C);
}


int main(int argc, char *argv[])
{
int rank, size;
int n1; // size of the square matrix
double *mat1;
double *matrix1;
double *matrix1_transposed;
double *result1_seq;
double *result1_parallel;
double start_time1, end_time1;
double start_time2, end_time2;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

n1 = pow(2,15);
mat1 = (double *)malloc(n1 * n1 * sizeof(double));
matrix1 = (double *)malloc(n1 * n1 * sizeof(double));
matrix1_transposed = (double *)malloc(n1 * n1 * sizeof(double));
result1_parallel = (double *)calloc(n1 * n1, sizeof(double));
result1_seq = (double *)calloc(n1 * n1, sizeof(double));

if (rank == 0)
{
    // Fill the matrix with random double precision floating point numbers between -1 and 1
    for (int i = 0; i < n1; i++)
    {
         for (int j = 0; j < n1; j++)
         {
              mat1[i * n1 + j] = (double)rand() / RAND_MAX * 2.0 - 1.0;
              matrix1[i * n1 + j] = mat1[i * n1 + j];
         }
    }
    transpose_matrix(n1, matrix1, matrix1_transposed);
    //print_matrix_to_file(mat1, n1, "a.txt"); //for debugging
    transpose_matrix(n1, matrix1, matrix1_transposed);
    //print_matrix_to_file(matrix1_transposed, n1, "b_transposed.txt"); //for debugging
    //print_matrix_to_file(matrix1, n1, "b.txt"); //for debugging
}

if(rank == 0){    
    start_time1 = MPI_Wtime();
    seq_mat_mult(mat1, matrix1, result1_seq, n1);
    end_time1 = MPI_Wtime();
    double time_spent1 = end_time1 - start_time1;
    printf("Time spent for sequential multiplying two %d by %d matrices is %f seconds. From rank: %d\n", n1, n1, time_spent1, rank);
}

MPI_Barrier(MPI_COMM_WORLD);

if(rank == 0){
    start_time2 = MPI_Wtime();
}
parallel_mat_mult(mat1, matrix1_transposed, result1_parallel, n1, size);
MPI_Barrier(MPI_COMM_WORLD);
if(rank == 0){
    end_time2 = MPI_Wtime();
    double time_spent2 = end_time2 - start_time2;
    printf("Time spent for parallel multiplying two %d by %d matrices is %f seconds. From rank: %d\n", n1, n1, time_spent2, rank);
}

MPI_Barrier(MPI_COMM_WORLD);
printf("Hello from %d of %d processors.\n", rank, size);
if(rank == 0){
    check_correctness(result1_seq, result1_parallel, n1);
}
MPI_Finalize();

// Free the memory
free(mat1);
free(matrix1);
free(matrix1_transposed);
free(result1_parallel);
free(result1_seq);

return 0;
}

