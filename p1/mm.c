/**
 * @file matrixmultiple.c
 *
 * @brief 
 *    Implementation of dot product of two nxn matrices.
 *
 * @author Athit Vue
 * @date 9/18/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Prototyping for multiply_matrices and print_matrix functions */
int *multiply_matrices (int *restrict first_matrix, int *restrict second_matrix, int matrix_size);
void print_matrix (int *matrix, int matrix_size);

/**
 * int main
 *
 * This is the main function. The multiplication of the matrices and
 * the printing of the resulting matrix has been decompose into its 
 * own functions. This function will ask the user to input R or I as 
 * the flag. If R is selected then resulting matrix is not printed. 
 * Else the resulting matrix will be printed.
 */
int main(void)
{
  // Declaration of variables
  char flag;
  int n;
  int *restrict matrix_A;
  int *restrict matrix_B;
  int *restrict matrix_C; // Result of matrix A and matrix B

  // Use scanf() to read from stdin
  scanf(" %c", &flag);
  scanf(" %d", &n);

  // Allocate memory for matrix A and matrix B using malloc
  matrix_A = malloc(sizeof(int) * n * n);
  matrix_B = malloc(sizeof(int) * n * n);

  // If flag is == I, then create matrix according to input.
  // Else, generate two random matrices...do not print resulting matrix.   
  if (flag == 'I')
  {
    // Use scanf to take in input for matrix A
    for (int row = 0; row < n; row++)
    {
      for (int col = 0; col < n; col++)
      {
        scanf(" %d", (matrix_A + row * n + col));
      }
    }
    // Use scanf to take in input for matrix B
    for (int row = 0; row < n; row++)
    {
      for (int col = 0; col < n; col++)
      {
        scanf(" %d", (matrix_B + row * n + col));
      }
    }

    // Need to get the resulting matrix (matrix C) and print it out
    matrix_C = multiply_matrices(matrix_A, matrix_B, n);
    print_matrix(matrix_C, n);
  }
  else
  {
    // Generate random numbers for the matrices seeding the generator
    // with the current time. Generate integers between 0 and 100.
    srand(time(0));
    int row, col;
    for (row = 0; row < n; row++)
    {
      for (col = 0; col < n; col++)
      {
        matrix_A[row*n + col] = rand() % 100;
        matrix_B[row*n + col] = rand() % 100;
      }
    }

    // Invoke multiply_matrices to multiply matrix A and matrix B
    matrix_C = multiply_matrices(matrix_A, matrix_B, n);
  }

  // Need to free all allocated memory from using malloc
  free(matrix_A);
  free(matrix_B);
  return 0;
}

/**
 * int multiply_matrices
 *
 * This function performs the dot product of two matrices - A and B. 
 *
 * @param first_matrix   The first matrix. 
 * @param second_matrix  The second matrix.
 * @return results       The resulting matrix.
 */
int *multiply_matrices(int *restrict first_matrix, int *restrict second_matrix, int matrix_size)
{
  int *restrict results = malloc(sizeof(int) * matrix_size * matrix_size);
  int i, j, k;
  //#pragma vector aligned
  for (i = 0; i < matrix_size; i++)
  {
    for (j = 0; j < matrix_size; j++)
    {
      for (k = 0; k < matrix_size; k++)
      {
        results[i*matrix_size + j] += first_matrix[i*matrix_size + k]
        * second_matrix[k*matrix_size + j];
      }
    }
  }
  return results;
}

/**
 * void print_matrix
 *
 * This function prints the matrix when the user
 * selects I as the flag from stdin.
 *
 * @param matrix      A pointer to the matrix to be printed.
 * @param matrix_size The size of the matrix. 
 */
void print_matrix(int *matrix, int matrix_size)
{
  for (int row = 0; row < matrix_size; row++)
  {
    //#pragma ivdep
    for (int col = 0; col < matrix_size; col++)
    {
      if (col != 0)
      {
        printf(" ");
      }
      printf("%d", matrix[row*matrix_size + col]);
    }
    printf("\n");
  }
}


