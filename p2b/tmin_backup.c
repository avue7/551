/** 
 * @file tmin.c
 *
 * @brief 
 *   This program performs a numerical integration using the 
 *   Trapezoidal Method. The upper and lower bound is accepted
 *   as the inputs and t to represent the minimum trapezoids
 *   need in order to be as accurate as we can to 14 significant
 *   figures of the true value of 4754.0192288588181366. The 
 *   function to be integrated is hardcoded below. This program 
 *   is a parallel program using MPI to calculate the definite
 *   integral. 
 *
 * @author Athit Vue
 * @date 10/01/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* Function Prototype */
long double estimate_integral(long double a, long double b, int n, 
  long double h);
long double abs_relative_true_error(long double approximate_value);
long double f(long double x);

/* For MPI */
void get_input(int my_rank, int comm_sz, long double* a,
  long double* b, int* n);
void create_mpi_type(long double* a, long double *b, int* n,
  MPI_Datatype* input_data_struct);

/* Global variables */
const long double TRUE_VALUE = 4754.0192288588181366;
const long double ACCEPTABLE_ERROR = .5E-14;

/**
 * int main
 *
 * This is the main function. This function accepts inputs as the 
 * left and right bound of the integral, and an input for the number
 * of trapezoids in order to search for the tmin. Once this function
 * is done calling the working functions, then it will print the 
 * integration result, the absolute relative true error, and the tmin
 * used. 
 */ 
int main(void)
{
  /* Declaration and initialization of variables */
  int comm_sz = 0; // Number of processes
  int my_rank = 0;
  int n = 0; // Trapezoids
  int local_n = 0;
  long double a = 0.0; // Left bound
  long double b = 0.0; // Right bound
  long double local_a = 0.0; // Local a
  long double local_b = 0.0; // Local b
  long double h = 0.0;
  long double local_approx_value = 0.0;
  long double total_approx_value = 0.0;
  double local_start_time = 0.0;
  double local_finish_time = 0.0;
  double local_elapsed_time = 0.0;
  double elapsed_time = 0.0;

  // Initialize MPI
  MPI_Init(NULL, NULL);
  // Get my rank
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // Get number of processes being used
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  // Get input
  get_input(my_rank, comm_sz, &a, &b, &n);
  // Block needed to start timer
  MPI_Barrier(MPI_COMM_WORLD);

  // Start the timer
  local_start_time = MPI_Wtime();

  /* Set the local variables for all ranks */
  h = (b-a)/n;
  local_n = n/comm_sz; // Same number of trapezoids for all
  local_a = a + my_rank*local_n*h;
  local_b = local_a + local_n*h;
  // Calculate the local_integral
  local_approx_value = estimate_integral(local_a, local_b, local_n, h);

  printf("LOCAL approx value = %.13e\n", (double) local_approx_value);
  // Add up all the local approximated values of each process
  MPI_Reduce(&local_approx_value, &total_approx_value, 1, MPI_DOUBLE,
    MPI_SUM, 0, MPI_COMM_WORLD);

  // Record the local finish time
  local_finish_time = MPI_Wtime();
  local_elapsed_time = local_finish_time - local_start_time;

  // Get the min local elapsed time and store it in variable elapsed_time
  MPI_Reduce(&local_elapsed_time, &elapsed_time, 1, MPI_DOUBLE, MPI_MIN,
    0, MPI_COMM_WORLD);
 
  /* PRINT OUT */
  if (my_rank == 0)
  {
    printf("Elapsed time = %e seconds\n", elapsed_time);
    printf("With n = %d trapezoids, our estimate\n", n); 
    printf("of the integral from %f to %f = %.13e\n",(double) a, (double) b,
      (double) total_approx_value);
    printf("true value = %.19e\n", (double) TRUE_VALUE);
  
    // Get the abs relative true error
    long double relative_true_error = 0.0;
    relative_true_error = abs_relative_true_error(total_approx_value);
    printf("absolute relative true error = %.19e\n", 
      (double) relative_true_error);
    if(relative_true_error > ACCEPTABLE_ERROR)
    {
      printf(" is NOT less than criteria = %.19e\n",
      (double) ACCEPTABLE_ERROR);
    }
  }
  // Clean up and shut down MPI
  MPI_Finalize();
}

/**
 * long double abs_relative_true_error
 * 
 * Calculates the absolute relative error.
 *
 * @param approximate_value The approximate_value attained from integration.
 * @return relative_true_error The relative true error.
 */
long double abs_relative_true_error(long double approximate_value)
{
  long double true_error = 0.0;
  true_error = TRUE_VALUE - approximate_value;
  long double relative_true_error = 0.0;
  relative_true_error = true_error / TRUE_VALUE;

  // Get the absolute value
  if (relative_true_error < 0)
  {
    relative_true_error = -relative_true_error;
  }

  return relative_true_error;
}

/**
 * long double estimate_integral
 * 
 * THis function uses the trapazoidal rule to approximate the result.
 * 
 * @param a The left bound of the integral.
 * @param b The right bound of the integral.
 * @t The number of trapezoids inputed.
 * @return approx The approximated result.
 */
long double estimate_integral(long double a, long double b, int n, 
  long double h)
{
  long double approx = 0.0;
  long double x_i = 0.0;
  int i = 0;

  approx = (f(a) + f(b))/2.0;
  for (i = 1; i <= n-1; i++)
  {
    x_i = a + i*h;
    approx += f(x_i);
  }
  approx = h*approx;
  return approx;
}

/**
 * long double f
 * 
 * This function calculates the integral for the given x
 * and returns it. 
 *
 * @param x The lower or upper bound to estimate the function.
 * @return Returns the estimated value.
 */
long double f(long double x)
{
  return -3*cosl(x/6) + sqrtl(powl((5*sinl(x/9))+(8*sinl(x/2)),4)) + 1;
}

/**
 * void get_input
 *
 * Gets the input from user and instead of broadcasting the three
 * different datatypes 3 times define a derived datatype and call 
 * MPI_Bcast only once.
 *
 * note: Borrowed from Pacheco's mpi_trap4.c example.
 *
 * @param a Pointer to the left endpoint
 * @param b Pointer to the right endpoint
 * @param n Pointer to the number of trapezoids
 */
void get_input(int my_rank, int comm_sz, long double* a,
  long double* b, int* n)
{
  MPI_Datatype input_data_struct;
  // Buld the mpi type
  create_mpi_type(a, b, n, &input_data_struct);
  
  if (my_rank == 0)
  {
    printf("Enter a, b, and n\n");
    // scan for user input for a, b, and t.
    scanf(" %LF %LF %d", a, b, n);
    // printf("Testing scan done!\n");
    printf("Running on %d cores.\n", comm_sz );
  }
  
  /* Now we can broadcast it only once */
  MPI_Bcast(a, 1, input_data_struct, 0, MPI_COMM_WORLD);
  // Free the data structure created
  MPI_Type_free(&input_data_struct);
}

/**
 * void create_mpi_type 
 *
 * Creates the data structure for the derived types.
 * 
 * note: Borrowed from Pacheco's mpi_trap4.c example.
 *
 * @param a Pointer to the left bound.
 * @param b Pointer to the right bound
 * @param n Pointer to the number of trapezoids.
 * @param input_data_struct Pointer to the data structure.
 */
void create_mpi_type(long double* a, long double* b, int* n,
  MPI_Datatype* input_data_struct)
{
  // Size of block
  int blocklengths_array[3] = {1,1,1};
  // Declare the MPI_types of each element in the array block
  MPI_Datatype types_array[3] = {MPI_LONG_DOUBLE,
    MPI_LONG_DOUBLE, MPI_INT};
  // Declare MPI_Aint type variable to store displacement
  MPI_Aint displacements_array[3] = {0};
  // Declare address variables to store the address
  MPI_Aint a_addr, b_addr, n_addr;

  /* Get the displacements for a, b, and n */
  MPI_Get_address(a, &a_addr);
  MPI_Get_address(b, &b_addr);
  MPI_Get_address(n, &n_addr);
 
  /* Set the address displacement for each type of datatypes*/
  displacements_array[1] = b_addr - a_addr;
  displacements_array[2] = n_addr - a_addr;
  
  // Create the derived data type
  MPI_Type_create_struct(3, blocklengths_array, displacements_array,
    types_array, input_data_struct);
  // Must commit to allow MPI to optimize the internal representation
  MPI_Type_commit(input_data_struct);
}
