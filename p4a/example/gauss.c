
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/resource.h>
#include <string.h>
#define DEBUG 1

void print_matrix(double * matrix, int n);
double d_abs(double d);
double partial_pivot(double * matrix, int n, int cur_col);

int main(int argc, char * argv[]) {
  if (argc != 2) {
    printf("wrong number of arguments\n");
    return 0;
  }

  int n = atoi(argv[1]);
  int n_plus_one = n + 1;

  //the augmented matrix (n X n+1)
  double * matrix = malloc(sizeof(double) * n * n_plus_one);
  //these will be reused for looping through the matrix
  int row, col;

  //take in the provided input for the augmented matrix
  if (n < 11) {
    for (row = 0; row < n; row++) {
      for (col = 0; col < n_plus_one; col++) {
        scanf(" %lf", (matrix + row * n_plus_one + col));
      }
    }
  //else n >= 5, so create a random matrix
  } else {
    srand48(time(0));

    double range = 2 * 1.0e6;
    //generate a matrix with numbers between -1.0e6 and 1.0e6
    for (row = 0; row < n; row++) {
      for (col = 0; col < n_plus_one; col++) {
        matrix[row * n_plus_one + col] = drand48() * range
            - (double) 1.0e6;
      }
    }
  }

  //store a backup of the matrix that won't be modified
  double * matrix_backup = malloc(n * n_plus_one * sizeof(double));
  memcpy(matrix_backup, matrix, n * n_plus_one * sizeof(double));

  //forward elimination:

  //this will be used for storing the x solutions
  double * x_array = malloc(sizeof(double) * n_plus_one);
  double scalar, diag;
  double * row_with_diag;

  int col_2;
  double * matrix_row;

  //for every column, get zeroes below the diagonal
  for (col = 0; col < n - 1; col++) {
    //do partial pivoting
    diag = partial_pivot(matrix, n, col);

    row_with_diag = matrix + col*n_plus_one;

    //go through the rows that are below the diagonal in this column
    for (row = col + 1; row < n; row++) {
      scalar = matrix[row*n_plus_one + col] / diag;

      matrix_row = matrix + row*n_plus_one;
      /*do the multiplication for the temp row and minus that from
      the row being annihilated so that we annihilate the row*/
      for (col_2 = col; col_2 < n_plus_one; col_2++) {
        matrix_row[col_2] -= scalar * row_with_diag[col_2];
      }
    }
  }

  //now done with forward elimination

  //back substitution:

  //start at the bottom row and go up
  for (row = n - 1; row >= 0; row--) {
    diag = matrix[row*n_plus_one + row];

    //store this solution in x_array
    x_array[row] = matrix[row*n_plus_one + n]/diag;

    int row_2;
    //for every row above the diagonal
    for (row_2 = 0; row_2 < row; row_2++) {
      //do: b -= x_just_found * matrix[row_2][column_of_diag]
      matrix[row_2*n_plus_one + n] -= x_array[row] *
          matrix[row_2*n_plus_one + row];
    }
  }

  //back substitution is now done, and result is in x_array

  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  printf("Cpu time: %f\n", usage.ru_utime.tv_sec + ((double)
      (usage.ru_utime.tv_usec / (double) 1000000)));

  printf("System time: %f\n", usage.ru_stime.tv_sec + ((double)
      (usage.ru_stime.tv_usec / (double) 1000000)));

  printf("maximum resident set size: %ld\n", usage.ru_maxrss);
  printf("soft page faults: %ld\n", usage.ru_minflt);
  printf("hard page faults: %ld\n", usage.ru_majflt);

  int i;
  if (n < 11 || DEBUG) {
    //print the original matrix
    print_matrix(matrix_backup, n);

    //print the solution as one row
    printf("Solution is: \n");
    for (i = 0; i < n; i++) {
      printf("%.10e", x_array[i]);

      if (i != n - 1)
        printf(" ");
    }
    printf("\n");
  }

  //residual vector:

  /*this will at first store the A*x, but will
  be reused to store the residual vector*/
  double * matrix_times_x = malloc(sizeof(double) * n);

  int k;
  //start doing the matrix multiplication to get A*x
  for (row = 0; row < n; row++) { //i
    double c_ij = 0;
    for (k = 0; k < n; k++) { //k
      c_ij += matrix_backup[row*n_plus_one + k] * x_array[k];
    }
    matrix_times_x[row] = c_ij;
  }

  double i_squared_norm = 0;
  /*store the residual vector in matrix_times_x by doing
  matrix_times_x - b, while simaltaneously
  calculating the i_squared_norm from this*/
  for (i = 0; i < n; i++) {
    matrix_times_x[i] -= matrix_backup[i*n_plus_one + n];
    i_squared_norm += matrix_times_x[i]*matrix_times_x[i];
  }

  i_squared_norm = sqrt(i_squared_norm);

  printf("i_squared norm: %.10e\n", i_squared_norm);

  //delete memory:

  free(matrix);
  free(matrix_backup);
  free(x_array);
  free(matrix_times_x);

  return 0;
}

/**
 * @brief   Prints a matrix.
 * @param   The matrix
 * @param   n, which is the size of the matrix (n X n+1)
 * @return  Nothing
 */
void print_matrix(double * matrix, int n) {
  int row, col;
  for (row = 0; row < n; row++) {
    for (col = 0; col < n + 1; col++) {
      if (col != 0)
        printf(" ");

      printf("%.10e", matrix[row * (n + 1) + col]);
    }
    printf("\n");
  }
}

/**
 * @brief   Takes the absolute value of a double.
 * @param   d, the double to take the absolute value of.
 * @return  The absolute value as a double.
 */
double d_abs(double d) {
  if (d < 0)
    return -d;
  else
    return d;
}

/**
 * @brief   Performs partial pivoting on the matrix given the column
 *            of the diagonal. The matrix provided will be modified.
 * @param   The matrix to perform on
 * @param   n, which is the size of the matrix (n X n+1)
 * @param   cur_col, the column that the diagonal will be on
 * @return  The best diagonal as a double.
 */
double partial_pivot(double * matrix, int n, int cur_col) {
  int row, n_plus_one = n + 1;
  double diag = matrix[cur_col*n_plus_one + cur_col];
  int best_row = cur_col;
  /*go through the rows that are at or below the  diagonal in
  this column looking for a better diagonal (partial pivoting)*/
  for (row = cur_col + 1; row < n; row++) {
    if (d_abs(matrix[row*n_plus_one + cur_col]) > d_abs(diag)) {
      best_row = row;
      diag = matrix[row*n_plus_one + cur_col];
    }
  }

  //if found another (better) row to swap with
  if (best_row != cur_col) {
    //swap the rows in place:

    int col;
    for (col = cur_col; col < n_plus_one; col++) {
      double tmp = matrix[best_row*n_plus_one + col];
      matrix[best_row*n_plus_one + col] = matrix[cur_col*n_plus_one + col];
      matrix[cur_col*n_plus_one + col] = tmp;
    }
  }

  return diag;
}
