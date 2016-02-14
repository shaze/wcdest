


#ifdef MPI
#include <mpi.h>
#else


#endif


// Constants used by calls by the MPI routines
#define ANSTAG 10
#define WORKTAG 11
#define MASTER 0




void transmitMPISlaveResponse();



// get the next block of work for MPI. NB: If MPI is not being used
// this still gets called, but it is trivial
// returns 0 if no more work to do, 1 if there is  work to do
// NB: Important side-effect -- sets global_i_beg, global_j_beg etc
// to indicate the bounds of the work to be done 
int mpi_get_task();



void do_MPImaster_cluster(WorkPtr work); 

void wcd_mpi_initialise(int argc, char *argv[]);

void wcd_mpi_cleanup();


