-- Options to import in the command line before execution. Need 
   to export the variables before you actually execute the program.

   * Block partitioning: (static scheduling)
     $ export OMP_SCHEDULE="static,5000"
     $ omp_sin_sum 2 10000

   * Cyclic partitioning: (static scheduling)
     $ export OMP_SCHEDULE="static,1"
     $ omp_sin_sum 2 10000

   * Guided partitioning: (guided scheduling)
     $ export OMP_SCHEDULE="guided"
     $ omp_sin_sum 2 10000
     
     Note: chunks are dynamically allocated to threads.
     In addition, size of the chunk is gradually decreased.
     Each allocation is approximately number of remaining
     iterations/number of threads. 

    
 **Additional Notes:
   * Static Schedule:
     - If chunksize is not specified, it defaults to number of 
       iterations/number of threads.
     - Chunks are statically assigned to threads (assignment is 
       done before the loop starts).

   * Dynamic Schedule: 
     - Done while the loop is executing.
     - When a thread finishes a chunk, it asks for another chunk. 
     - If chunksize is not specified, it defaults to 1.

   * Guided Schedule:
     - Chunks are dynamically allocated to threads.
     - In addition, size of the chunk is gradually decreased.
     - Each allocation is approximately number of remaining
       iterations/number of threads. 
   
   * Auto Schedule:
     - compiler or runtime system determine the schedule
   
   * Runtime Schedule:
     - an environment variable, OMP_SCHEDULE controls the schedule. 
     - Such as this program. 
