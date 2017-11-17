/* File:     omp_msglk.c
 *
 * Purpose:  Simulate message-passing using OpenMP.  This version uses
 *           an atomic directive and OpenMP locks to protect critical
 *           sections.
 *
 * Compile:  gcc -g -Wall -fopenmp -o omp_msglk omp_msglk.c queue_lk.c
 *           needs queue_lk.h
 * Usage:    ./omp_msglk <number of threads> <number of messages each 
 *                  thread sends>
 *
 * Input:    None
 * Output:   Source, destination and contents of each message received.
 *
 *
 * Notes:
 * 1.  DEBUG flag for more verbose output
 * 2.  This version uses locks to control access to the message queues.
 *
 * IPP:      Section 5.8.9 (pp. 248 and ff.)
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "queue_lk.h"
#include <time.h>
#define HIGHEST_PRIORITY 4

const int MAX_MSG = 10000;

void Usage(char* prog_name);
void Send_msg(struct queue_s* msg_queues[], int my_rank, 
      int thread_count, int msg_number, int priority);
void Try_receive(struct queue_s* q_p, int my_rank);
int Done(struct queue_s* q_p, int done_sending, int thread_count);

/*-------------------------------------------------------------------*/
int main(int argc, char* argv[]) {
   int thread_count;
   int send_max;
   struct queue_s** msg_queues;
   int done_sending = 0;

   if (argc != 3) Usage(argv[0]);
   thread_count = strtol(argv[1], NULL, 10);
   send_max = strtol(argv[2], NULL, 10);
   if (thread_count <= 0 || send_max < 0) Usage(argv[0]);

   msg_queues = malloc(thread_count*sizeof(struct queue_node_s*));
   memset(msg_queues, 0, thread_count*sizeof(struct queue_node_s*));

#  pragma omp parallel num_threads(thread_count) \
      default(none) shared(thread_count, send_max, msg_queues, done_sending)
   {
      int my_rank = omp_get_thread_num();
      int msg_number;
<<<<<<< HEAD
      srandom(my_rank + time(0));
=======
      srandom(time(0));
>>>>>>> 3cc9f4fae9cc6f8a2f95b4d8a534fa994256a5e8

      int priority = 0;

      msg_queues[my_rank] = Allocate_queue();

#     pragma omp barrier /* Don't let any threads send messages  */
                         /* until all queues are constructed     */

      for (msg_number = 0; msg_number < send_max; msg_number++) {
         priority = random() % 4 + 1; // Get random priority 1 - 4
         Send_msg(msg_queues, my_rank, thread_count, msg_number, priority);
         Try_receive(msg_queues[my_rank], my_rank);
      }
#     pragma omp atomic
      done_sending++;
#     ifdef DEBUG
      printf("Thread %d > done sending\n", my_rank);
#     endif

      while (!Done(msg_queues[my_rank], done_sending, thread_count))
      {
        //printf("In while loop for recive\n");
        Try_receive(msg_queues[my_rank], my_rank);
      }
      /* My queue is empty, and everyone is done sending             */
      /* So my queue won't be accessed again, and it's OK to free it */
      Free_queue(msg_queues[my_rank]);
      free(msg_queues[my_rank]);
   }  /* omp parallel */

   free(msg_queues);
   return 0;
}  /* main */

/*-------------------------------------------------------------------*/
void Usage(char *prog_name) {
   fprintf(stderr, "usage: %s <number of threads> <number of messages>\n",
         prog_name);
   fprintf(stderr, "   number of messages = number sent by each thread\n");
   exit(0);
}  /* Usage */

/*-------------------------------------------------------------------*/
void Send_msg(struct queue_s* msg_queues[], int my_rank, 
      int thread_count, int msg_number, int priority) {
// int mesg = random() % MAX_MSG;
   int mesg = msg_number;
   int dest = 0;

   // Make sure my_rank does not send to iteself
   while (1)
   {
     dest = random() % thread_count;
     if (dest != my_rank)
     {
       break;
     }
   }

   struct queue_s* q_p = msg_queues[dest];
   omp_set_lock(&q_p->lock);
   Enqueue(q_p, my_rank, mesg, priority);
   omp_unset_lock(&q_p->lock);
#  ifdef DEBUG
   printf("Thread %d > sent %d,%d to %d\n", my_rank, mesg, priority, dest);
#  endif
}  /* Send_msg */

/*-------------------------------------------------------------------*/
void Try_receive(struct queue_s* q_p, int my_rank) {
   int src, mesg, priority;
   int queue_size = q_p->enqueued - q_p->dequeued;
  // printf("queue_size in try_recieve %d\n", queue_size);

   if (queue_size == 0) return;
   else if (queue_size == 1) {
      omp_set_lock(&q_p->lock);
      Dequeue(q_p, &src, &mesg, &priority);  
      omp_unset_lock(&q_p->lock);
   } else{
      Dequeue(q_p, &src, &mesg, &priority);
   }
   printf("Thread %d > received %d,%d from %d\n", my_rank, mesg, priority, src);
}   /* Try_receive */

/*-------------------------------------------------------------------*/
int Done(struct queue_s* q_p, int done_sending, int thread_count) {
   int queue_size = q_p->enqueued - q_p->dequeued;
//   printf("Queue_size is %d\n", queue_size);
   if (queue_size == 0 && done_sending == thread_count)
      return 1;
   else 
      return 0;
}   /* Done */
