/**
 * @file omp_example
 *
 * @brief
 *       This program is an implementation of the trapezoid rule 
 *   using OpenMP. Users must specify the number of threads that 
 *   he or she wants from the command line argument. Then the 
 *   program will ask the user for the left bound, right bound, 
 *   and the number of trapezoids to calculate. OpenMP is then 
 *   implemented and the results will then be printed. 
 *
 * @author Athit Vue
 * @date 10/24/17
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/* Function Prototype */
long double estimate_integral(long double a, long double b, int n, int thread_count);
long double abs_relative_true_error(long double approximate_value);
long double f(long double x);
void print_my_rank();
int get_num_threads();

/* Global variables */
const long double TRUE_VALUE = 4754.0192288588181366L;
const long double ACCEPTABLE_ERROR = .5E-14L;

int main(int argc, char* argv[])
{
  /* Declaration and initialization of variables */
  int n = 0; // Trapezoids
  long double a = 0.0L;
  long double b = 0.0L;
  long double global_result = 0.0L;
  double start_time = 0.0;
  double finish_time = 0.0;
  double elapsed_time = 0.0;

  int thread_count;

  thread_count = strtol(argv[1], NULL, 10);
 
  printf("Enter a, b, and n\n");
  // scan for user input for a, b, and t.
  scanf(" %LF %LF %d", &a, &b, &n);

  start_time = omp_get_wtime();

  global_result = estimate_integral(a, b, n, thread_count);

  finish_time = omp_get_wtime();
  elapsed_time = finish_time - start_time;

  int check_thread_num = 0;
# pragma omp parallel num_threads(thread_count) 
  check_thread_num = get_num_threads();

  if (check_thread_num == 1)
  {
    printf("Running on %d thread.\n", check_thread_num);
  }
  else
  {
    printf("Running on %d threads.\n", check_thread_num);
  }
  printf("Elapsed time = %e seconds\n", elapsed_time);
  printf("With n = %d trapezoids, our estimate\n", n); 
  printf("of the integral from %f to %f = %.13e\n",(double) a, (double) b,
    (double) global_result);
  printf("true value = %.19e\n", (double) TRUE_VALUE);
  
  // Get the abs relative true error
  long double relative_true_error = 0.0;
  relative_true_error = abs_relative_true_error(global_result);
  printf("absolute relative true error = %.19e\n", 
    (double) relative_true_error);
  if(relative_true_error > ACCEPTABLE_ERROR)
  {
    printf(" is NOT less than criteria = %.19e\n",
    (double) ACCEPTABLE_ERROR);
  }
  return 0;
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
  long double true_error = 0.0L;
  true_error = TRUE_VALUE - approximate_value;
  long double relative_true_error = 0.0L;
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
  int thread_count)
{
  long double h = 0.0L;
  long double approx = 0.0L;
  long double x_i = 0.0L;
  int i = 0;

  h = (b-a)/(long double)n;
  approx = (f(a) + f(b))/2.0;

# pragma omp parallel for num_threads(thread_count) reduction(+: approx)
  for (i = 1; i <= n-1; i++)
  {
    x_i = a + i*h;
    approx += f(x_i);
  }
  approx = h*approx;

  /* Delete All within this block for better timing */
  printf("\n");
# pragma omp parallel num_threads(thread_count)
  print_my_rank();
  printf("\n");
  /*------------------------------------------------*/
  
  return approx;
}

/**
 * void print_my_rank
 *
 * Function to print what thread has finished working.
 *
 */
void print_my_rank()
{
  int my_thread_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();

  if (thread_count == 1)
  {
    printf("My rank %d of %d thread finished.\n", my_thread_rank, 
        thread_count);  
  }
  else
  {
    printf("My rank %d of %d threads finished.\n", my_thread_rank,
        thread_count); 
  }
}

/**
 * int get_num_threads()
 *
 * Function to get the number of threads that are working. 
 *
 * @return thread_count The number of threads.
 */
int get_num_threads()
{
  int thread_count = omp_get_num_threads();
  return thread_count;
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

