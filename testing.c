#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
  long double result;
  long double  x = 1;
  result = -3*cos(x/6) + sqrt(pow((5*sin(x/9))+(8*sin(x/2)),4)) + 1;
  printf ("results = %LF \n", result);
  long double long_double = -10.12345678912345678;
  if (long_double < 0)
  {
    long_double = -long_double;
  }

  long_double = floor(1000*long_double)/1000;

  printf ("absolute long double %.20e \n", (double)long_double);
}
