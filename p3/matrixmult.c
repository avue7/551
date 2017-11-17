/**
 * @file matrixmult.c
 * 
 * @brief
 *   Implementation of parallel matrix multiplication using MPI.
 *
 * @author Athit Vue
 * @date 10/27/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#define MAX_FORM_SIZE 4
#define DEBUG 0

/* Prototype of functions*/
void get_input(char* form, char* flag, int* n, int comm_sz);
void get_matrix_input(int* matrix_A, int* matrix_B, int n);
void randomize_matrices(int* matrix_A, int* matrix_B, int n);
void get_send_counts_displacements(int* send_counts, int* displacements,
  int comm_sz, int local_n, int n);
void print_matrix(int *resulting_matrix, int n);

/**
 * int main
 *
 * The main function. This function will ask the user for the form, flag,
 * size of matrix, A (nxn matrix), and B (nxn matrix). It will then use
 * MPI to perform the multiplication of the two matrices. If flag is set
 * to R then the resulting matrix will not be printed. Else if flag is set
 * to I, then the resulting matrix will be printed. The number of processors
 * running and elapsed time should be printed in both situation as well. 
 * 
 * @return 1 An error occurred.
 * @return 0 Program executed normally.
 */
int main()
{
  int my_rank = 0;
  int comm_sz = 0; // Number of processes
  int* matrix_A = NULL;
  int* matrix_B = NULL;
  int* matrix_C = NULL;
  char form[MAX_FORM_SIZE];
  char flag;
  int n = 0;
  int error = 0;
  double start_time = 0.0;
  double finish_time = 0.0; 
  double elapsed_time = 0.0;

  // Start the MPI
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  // Get my_rank
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Let rank 0 handle all inputs and randomization */
  if (my_rank == 0)
  {
    // Get input 
    get_input(form, &flag, &n, comm_sz);

    // Dynamically allocated space for matrix a and matrix b and c.
    // Only for process 0.
    matrix_A = malloc(sizeof(int)*n*n);
    matrix_B = malloc(sizeof(int)*n*n);
    matrix_C = malloc(sizeof(int)*n*n);

    // Zero all the matrices
    memset(matrix_A, 0, sizeof(int)*n*n);
    memset(matrix_B, 0, sizeof(int)*n*n);
    memset(matrix_C, 0, sizeof(int)*n*n);

    // If flag is I then we want to take in the matrix
    if (flag == 'I')
    {
      get_matrix_input(matrix_A, matrix_B, n);
      if (DEBUG) // Debugger
      {
        print_matrix(matrix_A, n);
        print_matrix(matrix_B, n);
      }
    }
    else if (flag == 'R') // Generate some random matrices
    {
      randomize_matrices(matrix_A, matrix_B, n);
      if (DEBUG) // Debugger
      {
        print_matrix(matrix_A, n);
        print_matrix(matrix_B, n);
      }
    }
    else // If flag != 'R' || flag != 'I', then invalid...exit.
    {
      printf("Invalid value for <flag> '%c', please enter 'R' or 'I'.\n",
        flag);
      error = 1;
      MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }
  
  // Check for error. If error, exit gracefully.
  if (error == 1)
  {
    if (my_rank == 0)
    {
      free(matrix_A);
      free(matrix_B);
      free(matrix_C);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);    
  }

  // Barrier to ensure all processes are here
  MPI_Barrier(MPI_COMM_WORLD);
  // Start the timer
  start_time = MPI_Wtime();

  /*######## INCLUDE ALL OVERHEADS STARTING HERE ###########*/
  
  // Broadcast the size of matrix from rank 0 to all others
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // Broadcast the form from rank 0 to all others 
  MPI_Bcast(&form, 4, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Get the local n size of matrix
  int local_n = 0;
  local_n = n / comm_sz;
  
  // Special case: when n cannot be divisible by comm_sz
  // (number of processors). The current processor should 
  // take in an extra row. Else, local_n is divisible by comm_sz
  int new_local_n = 0;
  if (n % comm_sz >= my_rank + 1)
  {
    new_local_n = local_n + 1;
  }
  else
  {
    new_local_n = local_n; 
  }

  // Declare each local matrices.
  int* local_matrix_A = NULL;
  int* local_matrix_B = NULL;
  int* local_matrix_C = NULL;

  // Will contain new_local_n rows after scatter call
  local_matrix_A = malloc(sizeof(int) * new_local_n * n);
  memset(local_matrix_A, 0, sizeof(int)*new_local_n*n);

  // Number of rows that will be send to all processes 
  int* send_counts = NULL;
  // Displacement of scatter
  int* displacements = NULL;
  
  if (my_rank == 0)
  {
    // Allocate space for send counts and displacements and memset to 0.
    send_counts = (int*) malloc(comm_sz*sizeof(int));
    displacements = (int*) malloc(comm_sz*sizeof(int));
    memset(send_counts, 0, sizeof(int)*comm_sz);
    memset(displacements, 0, sizeof(int)*comm_sz);

    // Get the send counts and displacements
    get_send_counts_displacements(send_counts, displacements, 
       comm_sz, local_n, n);
  }

  // Scatter the rows of matrix A to all processes. Use scatterv so that 
  // gaps can be allowed between messages. This will allow for irregular
  // message sizes and also so that data can be distributed to processes
  // in any order.
  MPI_Scatterv(matrix_A, send_counts, displacements, MPI_INT, local_matrix_A,
  new_local_n*n, MPI_INT, 0, MPI_COMM_WORLD);

  // Need to create memory for local matrix B in all worker procs
  local_matrix_B = malloc(sizeof(int)*n*n);
  local_matrix_C = malloc(sizeof(int)*new_local_n*n);
  memset(local_matrix_B, 0, sizeof(int)*n*n);
  memset(local_matrix_C, 0, sizeof(int)*new_local_n*n);

  if (my_rank > 0)
  {
    matrix_B = malloc(sizeof(int)*n*n);
    memset(matrix_B, 0, sizeof(int)*n*n);
  }
  // Bcast matrix B fully to worker processes
  MPI_Bcast(matrix_B, n*n, MPI_INT, 0, MPI_COMM_WORLD);
  local_matrix_B = matrix_B;

  // Do multiplication by form for each processes
  int i;
  int j;
  int k;
  // ijk form
  if (form[1] == 'j')
  {
    for (i = 0; i < new_local_n; i++)
    {
      for (j = 0; j < n; j++)
      {
        for (k = 0; k < n; k++)
        {
          local_matrix_C[(i*n)+j] += local_matrix_A[(i*n)+k]
          * local_matrix_B[(k*n)+j];
        }
      }
    }
  }
  // ikj
  else if (form[1] == 'k')
  {
    for (i = 0; i < new_local_n; i++)
    {
      for (k = 0; k < n; k++)
      { 
        for (j = 0; j < n; j++)
        {
          local_matrix_C[(i*n)+j] += local_matrix_A[(i*n)+k]
          * local_matrix_B[(k*n)+j];
        }
      }
    }
  }
  // kij
  else if (form[1] == 'i')
  {
    for (k = 0; k < n; k++)
    {
      for (i = 0; i < new_local_n; i++)
      {
        for (j = 0; j < n; j++)
        {
          local_matrix_C[(i*n)+j] += local_matrix_A[(i*n)+k]
          * local_matrix_B[(k*n)+j];
        }
      }
    }
  }
  else
  {
    printf("Somthing went wrong while trying to do the dot product!\n");
    if (my_rank == 0)
    {
      free(send_counts);
      free(displacements);
      free(matrix_A);
      free(matrix_C);
    }
    free(local_matrix_A);
    free(local_matrix_B);
    free(local_matrix_C);
    MPI_Abort(MPI_COMM_WORLD, 1);    
  }

  // Gather all the results from worker processes and add them to the 
  // result from process 0. Store the final result in process 0.
  // Use gatherV so that we can have irregular size. This is for when 
  // n is not divisible by number of processes.
  MPI_Gatherv(local_matrix_C, new_local_n*n, MPI_INT, matrix_C,
    send_counts, displacements, MPI_INT, 0, MPI_COMM_WORLD);
 
  if (my_rank == 0)
  {
    finish_time = MPI_Wtime();
    elapsed_time = finish_time - start_time;

    printf("elapsed time = %e seconds\n", elapsed_time);

    if (flag != 'R')
    {
      print_matrix(matrix_C, n);
    }
  }

  // Free all resources 
  if (my_rank == 0)
  {
    free(send_counts);
    free(displacements);
    free(matrix_A);
    free(matrix_C);
  }
  free(local_matrix_A);
  free(local_matrix_B);
  free(local_matrix_C);

  MPI_Finalize();

  return 0;
}

/**
 * void get_input
 *
 * This function is use to the form, flag, and n (size of matrix) from
 * the user. 
 */
void get_input(char* form, char* flag, int* n, int comm_sz)
{
  printf ("running on %d processors\n", comm_sz);
  scanf(" %s", form);
  scanf(" %c", flag);
  scanf(" %d", n);
}

/**
 * void get_matrix_input
 *
 * Takes in the inputs for matrix_A and matrix_B if 'I' was selected
 * as the flag.
 *
 * @param matrix_A Pointer to matrix A.
 * @param matrix_B Pointer to matrix B.
 * @param n The size of the matrix.
 */
void get_matrix_input(int* matrix_A, int* matrix_B, int n)
{
  int row;
  int col;
  // Take in inputs for matrix A
  for (row = 0; row < n; row++)
  {
    for (col = 0; col < n; col++)
    {
      scanf(" %d", (matrix_A + row * n + col));
    }
  }
  // Take in inputs for matrix B
  for (row = 0; row < n; row++)
  {
    for (col = 0; col < n; col++)
    {
      scanf(" %d", (matrix_B + row * n + col));
    }
  }   
}

/**
 * void randomize_matrices
 *
 * Function to randomize matrix A and matrix B.
 *
 * @param matrix_A Pointer to matrix A.
 * @param matrix_B Pointer to matrix B.
 * @param n The size of the matrix.
 */
void randomize_matrices(int* matrix_A, int* matrix_B, int n)
{
  srand(time(0)); // seed the random generator with current time
  int row;
  int col;
  for (row = 0; row < n; row++)
  {
    for (col = 0; col < n; col++)
    {
      matrix_A[(row*n) + col] = rand() % 100; // between 0 and 100
      matrix_B[(row*n) + col] = rand() % 100; // between 0 and 100
    }
  }
}

/**
 * void get_send_counts_displacements
 * 
 * Function to get the send_counts and displacement needed for 
 * call to scatter.
 * 
 * @param send_counts Pointer to send_counts.
 * @param displacements Pointer to displacements.
 * @param comm_sz The number of processors.
 * @param local_n The number of local n (partial size of matrix).
 * @param n The size of the matrix.
 */
void get_send_counts_displacements(int* send_counts, int* displacements,
  int comm_sz, int local_n, int n)
{
  int proc;
  // Get send counts for each processor (# of items to send).
  for (proc = 0; proc < comm_sz; proc++)
  {
    if (n % comm_sz >= proc+1)
    {
      send_counts[proc] = (local_n*n) + n;
    }
    else
    {
      send_counts[proc] = local_n*n;
    }
  }
  // Get displacements
  for (proc = 0; proc < comm_sz; proc++)
  {
    if (proc != 0)
    {
      displacements[proc] = displacements[proc-1] + send_counts[proc-1];
    }
    else
    {
      displacements[0] = 0;
    }
  }
}

/** 
 * void print_matrix
 *
 * Function to print the resulting matrix if "I" was inputted as
 * the flag by the user. 
 *
 * @param resulting_matrix Pointer to the the resulting matrix of
 *   the dot product matrix A and matrix B.
 * @param n The size of the matrix n X n.
 */
void print_matrix(int *resulting_matrix, int n)
{
  int row;
  int col;
  for (row = 0; row < n; row++)
  {
    for (col = 0; col < n; col++)
    {
      if (col != 0)
      {
        printf (" ");
      }
      printf("%d", resulting_matrix[(row*n) + col]);
    }
    printf("\n");
  }
}
