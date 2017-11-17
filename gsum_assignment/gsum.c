#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  int my_rank, comm_sz;
  float local_sum = 0;

  MPI_Init(NULL, NULL);

  if (strcmp(argv[0], "p") == 0)
  {
    /* Get # of pr0cesses */
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    // Get rank among all processes
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    float current_rank = my_rank;

    // Send rank to rank 0 if rank is not 0
    if (my_rank != 0)
    {
      MPI_Send(&current_rank, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
      for(int source = 1; source < comm_sz; source++)
      {
        MPI_Recv(&current_rank, 1, MPI_FLOAT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        local_sum += current_rank;
      }
      printf("The value of the global sum is: %f\n", local_sum);
    }
  }
  else if (strcmp(argv[0],"c") == 0)
  {
   
  }
  else 
  {
    printf("Please input p or c only!\n");
  }

  MPI_Finalize();
  return 0;
}
