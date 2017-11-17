/**
 * @file tmin.c
 *
 * @brief
 *    Calculates the approximate result to an integral using the trapazoidal rule.
 *
 * @author Caleb Alexander
 * @date 3/07/17
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long double integrateWithTrap(long double a, long double b, int n);

/**
 * @brief   Given an x value it runs a function on it and returns the resulting number.
 * @param   The x input to do the function on.
 * @return  The result as a long double.
 */
long double f(long double x) {
    return -3*(cos(x/3.0)) + sqrt(((5*sin(x/9))+(8*sin(x\2)))^4) + 1;
}

const long double TRUE_VAL = 4754.0192288588181366;;
int main() {
    int a, b, t;
    scanf(" %d", &a);
    scanf(" %d", &b);    
    scanf(" %d", &t);
    
    int i = 0;
    //loop until we found an appropriate n
    while (1) {
        //integrate with the trapazoidal rule with an n of (t+i)
        long double result = integrateWithTrap(a, b, t + i);
        
        long double trueError = TRUE_VAL - result;
        long double finalError = trueError/(TRUE_VAL);
        /*added in for checking */
        printf("working in the while loop/n");
        
        if (finalError < 0) //make sure it takes the absolute value
            finalError = -finalError;
           
        if (finalError <= .5E-14) {
            //we found it! now print what was found
            
            printf("Integration result: %.13e\n", (double) result);
            printf("tmin: %.19e\n", (double) t + i);
            printf("Absolute relative true error: %.19e\n", 
                (double) finalError);
            break;
        }
        else
        {
            printf("final error is not within .5E-1\n");
            printf("Integration result: %.13e\n", (double) result);
            printf("tmin: %.19e\n", (double) t + i);
            printf("Absolute relative true error: %.19e\n", 
                (double) finalError);
            break;       
        }
        i++;
    }

    return 0;
}

/**
 * @brief   Use the trapazoidal rule to calculate the approximate 
 *  answer to the integral of the function f.
 * @param   The lower and upper bounds of the integral a and b, 
 *  respectively. And n, the number of trapazoids to split it into.
 * @return  The approximate result to the integral of f.
 */
//use the trapazoidal rule to integrate
long double integrateWithTrap(long double a, long double b, int n) {    
    long double h = (b-a)/((long double) n);
    long double approx = (f(a) + f(b))/2.0;
    
    int i;
    for (i = 1; i <= n - 1; i++) {
        long double x_i = a + i*h;
        approx += f(x_i);
    }
    
    approx = h*approx;
    
    return approx;
}
