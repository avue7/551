/**
 * @file gauss-elim.c
 * 
 * @brief 
 *    This program performs the serial version of Gaussian Elimination
 *    with partial pivoting. 
 *
 * @author Athit Vue
 * @date 11/11/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>
#include <string.h>
#include <math.h>
#define DEBUG 0  // Set to 1 if you want to debug.

/* Prototyping */
void forward_elimination(double* augmented_matrix, int n);
void back_substitution(double* augmented_matrix, int n, double* solution_x); 
double partial_pivoting(double* augmented_matrix, int n, int col);
double abs_value(double value);
void do_matrix_multiplication(double* augmented_matrix_copy, 
  double* A_x, double* solution_x,  int n);
void get_L_sqaure_norm(double* augmented_matrix_copy, double* A_x,
  double* A_x_minus_b, double* L_square_norm, int n);
void print_matrix(double* augmented_matrix, int n);

/**
 * int main
 * 
 * This function is the main function. It expects one argument from the 
 * command line argument. This will be the size of the linear system. 
 * N will be nxn matrix A and n-vector b. The program will call 
 * functions foward_elimination() and back_substitution() to get the
 * solution x. After it has found the solution x, it will then print
 * the user cpu time, system cpu time, maximum resident set size, and 
 * minor and major page faults. If n < 11, then it will also print the 
 * original augmented matrix A before the forward elimination along with 
 * the solution x. An output of the L^2-norm will also be printed to check
 * the correctness of the solution x. 
 *
 * @param argc The command line argument count.
 * @param argv The string of characters of the command line argument.
 * @return 0 The program terminated gracefully.
 * @return 1 The program encounter an error.
 */
int main(int argc, char* argv[])
{
  int n; //Size of system.
  int row, col, i;
  double* augmented_matrix;
  double* augmented_matrix_copy; // Need if n > 11
  double* solution_x; // Need to store x solutions in an array
  double* A_x; // For storing A*x 
  double* A_x_minus_b; // For storing Ax-b
  double L_square_norm;

  // Error if more than one argument is given
  if (argc != 2)
  {
    if (argc < 2)
    {
      printf("Error, one argument is expected.\n");
    }
    else
    {
      printf("Error, only one argument is expected.\n");
    }
    return 1;
  }

  // Convert the first command line argument to an int.
  n = atoi(argv[1]); 
  
  /* Allocating augmented matrix; nxn matrix A and n-vector b.
     Matrix picks up and extra col for augmented matrix. 
     Example: n = 2: 1 2 | 5
                     3 4 | 6
  */ 
  augmented_matrix = malloc(sizeof(double) * n * (n + 1));

  // If n < 11, then read the augmented matrix for A|b from 
  // standard input, one row per line. 
  if (n < 11) 
  {
    for (row = 0; row < n; row++)
    {
      for (col = 0; col < (n + 1); col++)
      {
        if (scanf(" %lf", (augmented_matrix + row * (n + 1) + col)) == 1)
        {
          continue;
        }
      }
    }    
  }
  // N > 10, fill A|b with random numbers using drand48 random number 
  // generator.
  else
  {
    // Scale results to range = [-1.0e6, 1.0e6).
    double result_range = 2 * 1.0e6;
    // Seed the generator once with current time.
    srand48(time(0));
    // Fill A|b with random numbers.
    for (row = 0; row < n; row++)
    {
      for (col = 0; col < (n+1); col++)
      { 
        augmented_matrix[row * (n + 1) + col] = drand48()
        * result_range - (double) 1.0e6;
      }
    }    
  }
  
  // We will be needing the original augmented matrix for later
  // use so store a copy of it.
  augmented_matrix_copy = malloc(sizeof(double) * n * (n + 1));
  memcpy(augmented_matrix_copy, augmented_matrix, 
    sizeof(double) * n * (n+1));   
  
  // Do forward elimination on augmented matrix.
  forward_elimination(augmented_matrix, n);

  // Do back substitution to get the results.
  solution_x = malloc(sizeof(double)*(n+1));
  back_substitution(augmented_matrix, n, solution_x);

  // Need to make a call to getrusage once we have the solution x
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  printf("\n");
  printf("User CPU time: %f\n", usage.ru_utime.tv_sec + ((double)
    (usage.ru_utime.tv_usec / (double) 1000000)));
  printf("System CPU time: %f\n", usage.ru_stime.tv_sec + ((double)
    (usage.ru_utime.tv_usec / (double) 1000000)));
  printf("maximum resident set size: %ld\n", usage.ru_maxrss);
  printf("soft page faults: %ld\n", usage.ru_minflt);
  printf("hard page faults: %ld\n", usage.ru_majflt);
  printf("\n");

  // If n > 11, print the augmented matrix and solution x
  if (n < 11 || DEBUG)
  {
    printf("The original augmented matrix is:\n");
    print_matrix(augmented_matrix_copy, n);
    printf("\n");
    // Print the solution x
    printf("The solution x is:\n");

    for (i = 0; i < n; i++)
    {
      printf("%.10e", solution_x[i]);
      if (i != n-1)
      {
        printf(" ");
      }
    }
    printf("\n");
    printf("\n");
  }

  // Compute the square root of the sum of the squares of the 
  // residual vector Ax-b:

  // Need to allocate memory for Ax and Ax-b.
  A_x = malloc(sizeof(double)*n);
  A_x_minus_b = malloc(sizeof(double)*n);
  
  // Do matrix multiplication for A*x
  do_matrix_multiplication(augmented_matrix_copy, A_x, solution_x, n);

  // Get the L^2 norm
  get_L_sqaure_norm(augmented_matrix_copy, A_x, A_x_minus_b,
    &L_square_norm, n);

  // Print the L square norm
  printf("L^2-norm: %.10e\n", L_square_norm);

  // Free all resources 
  free(augmented_matrix);
  free(augmented_matrix_copy);
  free(solution_x);
  free(A_x);
  free(A_x_minus_b);
 
  return 0;  
}

/**
 * void forward_elimination
 *
 * This function performs forward elimination with partial pivoting.
 * Before every forward elimination process, the pivot (diagonal) will
 * be checked to see if it is the absolute largest of the values between
 * it and the lower values in the column. If one of the lower values is 
 * absolutely larger than the current pivot, then the their rows will be 
 * swapped when the call to partial_piviot() returns. 
 * 
 * @param augmented_matrix Pointer to the augmented matrix.
 * @param n The size of the augmented matrix (nx(n+1)).
 */
void forward_elimination(double* augmented_matrix, int n)
{
  int col, row_below, col_below;
  double factor; // Factor used for annihilation.
  double diag; // The diagonol of the matrix; the pivot term. 
  double* diag_row; // Pointer to the row with the diagonol.
  double* aug_matrix_row; // Pointer to augmented matrix row

  for (col = 0; col < n-1; col++)
  {
    // Before every forward elimination process, find the diagonal
    // of the best row and perform swap if necessary at function. 
    diag = partial_pivoting(augmented_matrix, n, col);
    // Set row to perform forward elimination
    diag_row = augmented_matrix + col * (n+1);
    
    // For rows below the diag in the column do annihilation.
    for (row_below = col + 1; row_below < n; row_below++)
    {
      factor = augmented_matrix[row_below*(n+1)+col] / diag;
      aug_matrix_row = augmented_matrix + row_below*(n+1);

      for (col_below = col; col_below < (n+1); col_below++)
      {
        aug_matrix_row[col_below] -= factor*diag_row[col_below];
      }
    }
  }     
}
 
/**
 * void back_substitution
 *
 * This function performs a back substitution on the augmented matrix
 * after the forward elimination processes has all completed. 
 * 
 * @param augmented_matrix Pointer to the augmented_matrix.
 * @param n The size of the system (nx(n+1)).
 * @param solution_x Pointer to the solution x.
 */
void back_substitution(double* augmented_matrix, int n, double* solution_x)
{
  int row_1, row_2;
  double diag;
  // Start from bottom to top
  for (row_1 = n-1; row_1 >= 0; row_1--)
  {
    // Get the diagonal
    diag = augmented_matrix[row_1*(n+1)+row_1];

    // Store the solution x 
    solution_x[row_1] = augmented_matrix[row_1*(n+1)+n] / diag;

    for (row_2 = 0; row_2 < row_1; row_2++)
    {
      augmented_matrix[row_2*(n+1)+n] -= solution_x[row_1]
      * augmented_matrix[row_2*(n+1)+row_1];
    }
  }
}

/**
 * double partial_pivoting
 *
 * This function is used to return the diagonol used for 
 * finding the solution of Gaussian Elimination. Partial pivoting
 * is used to avoid division by zero and reduces round off error. 
 * Before each forward elimination step, the best row from the 
 * diagonal down is selected to be used for the next elimination 
 * step. The best row is the one with the largest absolute value 
 * in the column being worked on. This function will return the 
 * best diagonal of the best row.  
 * 
 * @param augmented_matrix Pointer to the augmented matrix.
 * @param n The size of the augmented matrix.
 * @param col The current col we are working on.
 * @return diag The diag of the best row we have found. 
 */
double partial_pivoting(double* augmented_matrix, int n, int col)
{
  int row;
  int best; // Best row (largest absolute value)
  double diag;

  best = col;
  diag = augmented_matrix[col * (n+1) + col];
  
  // For each column, look for the largest absolute value at and 
  // below the diagonal. 
  for (row = col + 1; row < n; row++)
  {
    if (abs_value(augmented_matrix[row*(n+1)+col]) > abs(diag))
    {
      best = row;
      diag = augmented_matrix[row*(n+1)+col];
    }
  }

  // If the current col is not the best row, then find the the best 
  // row and swap it. 
  if (best != col)
  {
    int new_col;
    double temp;
    for (new_col = col; new_col < n+1; new_col++)
    {
      temp = augmented_matrix[best*(n+1)+new_col];
      augmented_matrix[best*(n+1)+new_col] = 
      augmented_matrix[col*(n+1)+new_col];
      augmented_matrix[col*(n+1)+new_col] = temp;
    }
  }
  return diag;
}

/**
 * double abs_value
 *
 * This function is used to return the absoulute value of a double 
 * value that is passed to it. This function is a helper function 
 * for the partial_pivoting function. 
 * 
 * @param value The value to check for asoluteness.
 * @return abs_value The absolute value of the double passed in. 
 */
double abs_value(double value)
{
  double abs_value = 0;
  if (value < 0)
  {
    abs_value = -value;
  }
  else
  {
    abs_value = value;
  }

  return abs_value;
}

/** 
 * void do_matrix_multiplication
 * 
 * This function performs the matrix multiplication to get
 * Ax. This will be used to help with computing Ax-b for computing
 * the residual vector. 
 *
 * @param augmented_matrix_copy Pointer to the original augmented matrix.
 * @param A_x Pointer to the array A_x.
 * @param solution_x Pointer to the solution x.
 * @param n The size of the system (nx(n+1)).
 */
void do_matrix_multiplication(double* augmented_matrix_copy, 
  double* A_x, double* solution_x,  int n)
{
  int i, k;
  for (i = 0; i < n; i++)
  {
    double sum = 0;
    for (k = 0; k < n; k++)
    {
      sum += augmented_matrix_copy[i*(n+1)+k] * solution_x[k];
    }
    A_x[i] = sum;
  }
}

/**
 * void get_L_square_norm
 * 
 * This function computes the L square norm. This L square norm
 * will be used by the program to check for the correctness of the
 * results (solution x) that is computed. 
 * 
 * @param augmented_matrix_copy Pointer to the original augmented matrix.
 * @param A_x Pointer to the array A_x (A times x).
 * @param L_square_norm Pointer to the address of L_square_norm.
 * @param n The size of the system (nx(n+1)).
 */
void get_L_sqaure_norm(double* augmented_matrix_copy, double* A_x,
  double* A_x_minus_b, double* L_square_norm, int n)
{
  int row;
  for (row = 0; row < n; row++)
  {
    A_x_minus_b[row] = A_x[row] - augmented_matrix_copy[row*(n+1)+n];
    *L_square_norm += A_x_minus_b[row] * A_x_minus_b[row];
  }
  
  *L_square_norm = sqrt(*L_square_norm);
}

/**
 * void print_matrix
 *
 * Prints the augmented matrix for A|b, one row per line.
 *
 * @param augmented_matrix The augmented matrix
 * @param n The size of the matrix system (nx(n+1)).
 */
void print_matrix(double* augmented_matrix, int n)
{
  int row, col;
  for (row = 0; row < n; row++)
  {
    for (col = 0; col < (n+1); col++)
    {
      if (col != 0)
      {
        printf(" ");
      }
      printf("%.10e", augmented_matrix[row*(n+1)+col]);
    }
    printf("\n");   
  }
}
