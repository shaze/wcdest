

#include "common.h"
#include "wcd.h"
#include "pthread.h"
#include "mpistuff.h"

extern int num_threads;
extern int global_i_beg, global_i_end, global_j_beg, global_j_end;
extern SeqPtr seq;




#ifdef PTHREADS
pthread_mutex_t get_task_mutex;
pthread_mutex_t write_union_mutex;
pthread_mutex_t invert_mutex;
pthread_mutex_t dump_mutex;
pthread_mutex_t find_parent_mutex;


int parallel_get_task(WorkPtr work) {

  int range, size, extra;
  int cursize = 0, diagonal;
  int d; // number of ways to divide i range
  static int i_cur=MAX_SEQ_LEN; // any arb very big number would do
 
  // This is set to false by a thread which finds local work
  // exhausted + no more remote work

  if (num_threads==1) 
    d=1;
  else
    d=8*num_threads;

  //DBg  printf("T%d i_cur=%d, d=%d\n",work->thread_num,i_cur,d);
  pthread_mutex_lock(&get_task_mutex);
  if (i_cur >= d) {  /* local tasks complete */
    if (mpi_get_task(work)) {// get remote work
      //mpi_get_task sets the global_X_Y variables
      i_cur =0;
    } else {   // no more remote work
      pthread_mutex_unlock(&get_task_mutex);
      i_cur=MAX_SEQ_LEN;
      return 0;
    }
  }
  range = global_i_end - global_i_beg;
  size = range / d;
  extra = range % d;

  // first "extra" blocks get one more row

  if (i_cur < extra) {
    work->i_beg = global_i_beg+i_cur * (size+1);
    work->i_end = work->i_beg + size + 1;
  } else {
    work->i_beg = global_i_beg+i_cur * size + extra;
    work->i_end = work->i_beg + size;
  }

  work->j_beg = global_j_beg;
  work->j_end = global_j_end;

  i_cur++;
  pthread_mutex_unlock(&get_task_mutex);
  return 1;
}
#endif

#ifdef PXXXXTHREADS
int parallel_get_task(WorkPtr work) {

  // This code originally written by Richard Starfield.
  // work   is the work block that this *thread* must do
  // -- this is the output from the function
  //    All the work to be done bythe process is divided into blocks
  //    (a 2-D matrix)
  //    The block (i_cur, j_cur) is returned
  //  The idea was that by allocating work in blocks that did not share
  //  sequences there would be less memory contention. However, experiments
  //  showed it didn't make much difference and once we introduced MPI+
  //  pthreads sharig it would have to become more complex. I've left this
  //  in incase anyone wants to experiment further
  


  int cursize = 0;
  int d; // number of ways to divide i and j ranges
  static int i_cur=MAX_SEQ_LEN;// i.e. some very big number
  static int j_cur=0;
  static int diagonal=MAX_SEQ_LEN;  // i.e. some very big number
  static int range, size, extra;


  if (num_threads==1) 
    d=1;
  else
    d=8*num_threads;
  // issue work by setting pointers


  pthread_mutex_lock(&get_task_mutex);

  if (diagonal >= d) {  /* local tasks complete */
    //printf("Asking for new owrk\n");
    if (mpi_get_task(work)) {
      range = global_i_end - global_i_beg;
      size = range / d;
      extra = range % d;
      i_cur = j_cur = 0;
      diagonal=0;
      //printf("Got new work: range=%d, size=%d; i=(%d,%d)\n",range,size,global_i_beg,global_i_end); 
    } else {
      pthread_mutex_unlock(&get_task_mutex);
      i_cur=MAX_SEQ_LEN;
      return 0;
    }
  }
  printf("T%d-%d: i_cur=%d, diagonal=%d, d=%d\n",myid,work->thread_num,i_cur,diagonal,d);


  // first "extra" blocks get one more row
  if (i_cur < extra) {
    work->i_beg = global_i_beg+i_cur * (size+1);
    work->i_end = work->i_beg + size + 1;
  } else {
    work->i_beg = global_i_beg+i_cur * size + extra;
    work->i_end = work->i_beg + size;
  }

  if (j_cur < extra) {
    work->j_beg = global_j_beg+j_cur * (size+1);
    work->j_end = work->j_beg + size + 1;
  } else {
    work->j_beg = global_j_beg+j_cur * size + extra;
    work->j_end = work->j_beg + size;
  }
  //DBg   
  printf("T%d-%d;   i %d  %d;    j %d %d\n",myid,work->thread_num,work->i_beg,work->i_end,work->j_beg,work->j_end);
  // Modified by SH  
  i_cur++;
  j_cur++;
  if (j_cur >= d) {
    j_cur = diagonal+1;
    diagonal++;
    i_cur = 0;
  }
  pthread_mutex_unlock(&get_task_mutex);
  return 1;
}
#endif

#if PTHREADS
// This is the main routine when  parallel_do_cluster is called
void parallel_do_cluster_servant(WorkPtr work) {
  int i_cur, j_cur;

  while (parallel_get_task(work)) {
    //DBg    printf("About to pairwise cluster\n");
    do_pairwise_cluster(work);
  }
  //printf("Thread %d finished\n",work->thread_num);
}
#endif


