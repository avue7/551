/**
 * Name:        Athit Vue
 * Program:     treesum.c
 * Description: A simple program outputs the treesum in a list of parallel 
 *              communication.
 **/ 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char* argv[])
{
  int rank;
  int divisor = 2;
  int rank_difference = 1;  
  int MAX_PROCS = atoi(argv[1]);
  int* procs_array[MAX_PROCS];

  /**************OPTIONAL FOR THIS PROGRAM I THINK***************/
  
  //initialize (random if want)  numbers for procs
  for (int i = 0; i < MAX_PROCS; i++)
  {
    int data = rand() % 15;
    procs_array[i] = data;
    //printf("processor %i data stored = %d\n", i, active_procs[i]);
  }
  /**************************************************************/

  /* We are assuming cores are powers of 2. This for loop will give us the 
   * height of the tree, the timestamp. Inside this loop we will perform
   * a series of checks for the sender and the reciever. It is important
   * to note here that I have implemented to make the value stored at
   * processor rank to be NULL, if it is a sender. This will allow me to 
   * check for a processor that has finished it's task in computing the sum
   * that it was responsible for doing. 
   */
  for (int i = 0; i < log2(MAX_PROCS); i++)
  {
    printf("timestamp %d: \n", i);
    for (int rank = 0; rank < MAX_PROCS; rank++) 
    { 
      if (procs_array[rank] !=  NULL)
      {
        // printf("---rank active =  %d--- \n", rank);
        /* check for the reciever, this processor still has task to do */
        if (rank % divisor == 0)
        {
          printf("%d recieves from %d\n", rank, rank+rank_difference);
        }
        /* check for the sender, this processor is done in the current
           timestamp. */
        else
        {
          printf("%d sends to %d\n", rank, rank-rank_difference);
          procs_array[rank] = NULL;
        }
      }
    } 
    divisor *= 2;
    rank_difference *= 2;
  }
}
