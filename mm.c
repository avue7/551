/**
 * @file mm.c
 *
 * @brief
 *    Multiplies two n X n matrices.
 *
 * @author Caleb Alexander
 * @date 2/21/17
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printMatrix(int * matrix, int n);
int * multiplyMatrix(int * matrix1, int * matrix2, int n);

int main() {
  char flag;
  scanf(" %c", &flag);
  int n; //each matrix is n X n
  scanf(" %d", &n);

  int * matrix1;
  int * matrix2;

  matrix1 = malloc(sizeof(int) * n * n);
  matrix2 = malloc(sizeof(int) * n * n);

  int * matrix_result;
  //generate random matrices
  if (flag == 'R') {
    srand(time(0));
    int row, col;
    //generate the 2 matrices with numbers between 0 and 100
    for (row = 0; row < n; row++) {
      for (col = 0; col < n; col++) {
        matrix1[row * n + col] = rand() % 100;
        matrix2[row * n + col] = rand() % 100;
      }
    }
    //multiply the matrices
    matrix_result = multiplyMatrix(matrix1, matrix2, n);
  //else take in matrix input
  } else {
    int row, col;
    //take in the input for matrix one
    for (row = 0; row < n; row++) {
      for (col = 0; col < n; col++) {
        scanf(" %d", (matrix1 + row * n + col));
      }
    }

    //take in the input for matrix two
    for (row = 0; row < n; row++) {
      for (col = 0; col < n; col++) {
        scanf(" %d", (matrix2 + row * n + col));
      }
    }
    //multiply the matrices
    matrix_result = multiplyMatrix(matrix1, matrix2, n);
    printMatrix(matrix_result, n);
  }

  //release all allocated memory:
  free(matrix1);
  free(matrix2);
  free(matrix_result);
  return 0;
}

/**
 * @brief   Prints a matrix.
 * @param   The matrix and n which is the size of the matrix (n X n).
 * @return  Nothing
 */
void printMatrix(int * matrix, int n) {
  int row, col;
  for (row = 0; row < n; row++) {
    for (col = 0; col < n; col++) {
      if (col != 0)
        printf(" ");

      printf("%d", matrix[row*n + col]);
    }
    printf("\n");
  }
}

/**
 * @brief   Multiplies 2 matrices and returns the resulting matrix.
 * @param   The two matrices to multiply together, and n which
 *            is the size of the matrix (n X n).
 * @return  The resulting matrix.
 */
int * multiplyMatrix(int * matrix1, int * matrix2, int n) {
  int * matrix_result = malloc(sizeof(int) * n * n);
  int row, col, i;
  //start doing the matrix multiplication
  for (row = 0; row < n; row++) {
    for (col = 0; col < n; col++) {
      for (i = 0; i < n; i++) {
        matrix_result[row*n + col] += matrix1[row*n + i] *
            matrix2[i*n + col];
      }
    }
  }
  return matrix_result;
}
