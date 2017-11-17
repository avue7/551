/** 
 * @file tmin.c
 *
 * @brief 
 *   This program performs a numerical integration using the 
 *   Trapezoidal Method. The upper and lower bound is accepted
 *   as the inputs and t to represent the minimum trapezoids
 *   need in order to be as accurate as we can to 14 significant
 *   figures of the true value of 4754.0192288588181366. The 
 *   function to be integrated is hardcoded below. 
 *
 * @author Athit Vue
 * @date 10/01/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Function Prototype */
long double estimate_integral(long double a, long double b, int t);
long double abs_relative_true_error(long double approximate_value);
long double f(long double x);

/* Global variables */
const long double true_value = 4754.0192288588181366L;
const long double acceptable_error = .5E-14L;

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
  long double a = 0.0L;
  long double b = 0.0L;
  // scan for user input for a, b, and t.
  scanf(" %LF", &a);
  scanf(" %LF", &b);
  int t = 0;;
  scanf(" %d", &t);
  // printf("Testing scan done!\n");
  
  int iter = 0; // Using this as iteration
  
  // Printing a new line for easier readability
  printf("\n");
  printf("################# RESULTS #####################\n");
  // Using infinite while loop
  while (1)
  {
    long double approximate_value = 0.0L;
    long double relative_true_error = 0.0L;

    approximate_value = estimate_integral(a, b, t + iter);

    // printf("returned successfully from estimate_integral \n");
 
    // Get the abs relative true error
    relative_true_error = abs_relative_true_error(approximate_value);
    
    if (relative_true_error <= acceptable_error)
    {
      printf("SUCCESS, found a possible tmin!\n");
      printf("1. Result from integration: %.13e\n", 
        (double) approximate_value);
      printf("2. Absolute relative true error: %.19e\n",
        (double) relative_true_error);
      printf("3. Tmin is: %.19e\n", (double) t + iter);
      break;
    }
    /*else
    {
      printf("FAILED, could not find tmin within acceptable error!\n");
      printf("1. Result from integration: %.13e\n", (double) approximate_value);
      printf("2. Absolute relative true error: %.19e\n",
        (double) relative_true_error);
      printf("3. Tmin is: %.19e\n", (double) t + iter);
      break; 
    }*/ 
    iter += 100;          
  }  
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
  true_error = true_value - approximate_value;
  long double relative_true_error = 0.0;
  relative_true_error = true_error / true_value;

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
long double estimate_integral(long double a, long double b, int t)
{
  long double h = 0.0;
  long double approx = 0.0;
  h = (b-a) / ((long double)t);
  approx = (f(a) + f(b))/2.0;

  for (int i = 1; i <= t-1; i++)
  {
    long double x_i = a + i*h;
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
