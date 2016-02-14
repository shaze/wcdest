


#include "common.h"
#include "wcd.h"
#ifdef PTHREADS
#include "pthread.h"
#include "math.h"
#include "pthreadstuff.h"
#endif



#ifdef MPIORIENT
#define MPIBUFRECSZ 4
#else
#define MPIBUFRECSZ 2
#endif



extern int myid, numprocs;



#ifdef MPI

#include <strings.h>
#include "mpi.h"

#define ANSTAG 10
#define WORKTAG 11
#define MASTER 0

int    merr;
static char error_string[MPI_MAX_ERROR_STRING];
extern int global_i_beg, global_i_end, global_j_beg, global_j_end;
extern int num_seqs, data_size;
extern SeqPtr *seq;
extern SeqPtr data;
extern SeqInfoPtr   seqInfo;
extern UnionFindPtr tree;
extern ProgOptionsType  prog_opts;

int32_t *buffer=NULL;


void mpierr(int error_code) {
  return;
  if (error_code == MPI_SUCCESS) return;
  int length_of_error_string;
  MPI_Error_string(error_code, error_string, &length_of_error_string);
  fprintf(stderr, "MPI process problem: %3d: %s\n", myid, error_string);
}



//---------- parallelisation


void transmitMPISlaveResponse() {
#ifdef MPI

  int      i, j,curr;
  curr=0;
  // Extract out the non-trivial clusters
  // printf("Slave %d responding\n",myid);
  for(i=0; i<num_seqs; i++) {
    if (tree[i].cluster == i && tree[i].next != -1)
      for(j=tree[i].next; j != -1; j=tree[j].next) {
	buffer[curr]=j;                 // seq id
	curr++;
	buffer[curr]=i;   // cluster
	curr++;
#ifdef MPIORIENT
	buffer[curr]=tree[j].match;
	curr++;
	buffer[curr]=tree[j].orient;
	curr++;
#endif
      }
  }
  assert(curr<MPIBUFRECSZ*num_seqs);
  MPI_Send(buffer, curr, MPI_INT, MASTER, ANSTAG, MPI_COMM_WORLD); 
#endif  
}




//---------- parallelisation

#ifdef MPI
void MergeSlaveClusterTable(int32_t *buffer, int rlen) {
  // Get the slave's cluster table and add to ours
  int i, id,  rep, rep_root, last;
  int q;
  assert(rlen%MPIBUFRECSZ==0);
  assert(rlen<MPIBUFRECSZ*num_seqs);
  last = -1;
  for(i=0; i<rlen; ) {
    id = buffer[i];
    i++;
    rep= buffer[i];
    i++;
#ifdef MPIORIENT
    int witness, orient;
    witness = buffer[i];
    i++;
    orient = buffer[i];
    i++;
#endif
    PLOCK(&write_union_mutex);
    PLOCK(&find_parent_mutex);
    if (last != rep) {
      last = rep;
      rep_root = find_parent(rep);
    }
    id = find_parent(id);
    if (id == rep_root) {
      PUNLOCK(&find_parent_mutex);
      PUNLOCK(&write_union_mutex);
      continue;
    }
    // if (orient) invert_orientation(j);
    q = tree[rep_root].last;
    tree[q].next = id;
    tree[rep_root].last = tree[id].last;
    tree[id].cluster = rep_root;
    if (IS_FLAG(FIX,id)) SET_FLAG(FIX,rep_root);
    assert(tree[tree[rep_root].last].next<0);
    PUNLOCK(&write_union_mutex);
    PUNLOCK(&find_parent_mutex);
  }
}


#endif



                 

void recalibrateSeqStruct() {
  int i;
  SeqPtr curr;

  curr = data;
  for(i=0; i<num_seqs; i++) {
    seq[i] = curr;
    curr = seq[i]+b2e(seqInfo[i].len);
  }
  assert(seq[0] == data);
}
#endif




void handleMPISlaveSetup(int *num_seqs) {
#ifdef MPI
  int i,j,k,nlen,shift,maxlen;
  double total_load=0.0,ave,curr;
  long tot_size=0, tot;
  int metadata[2];

  //printf("Slave %d setting up \n",myid);

  // Find out how many sequences there are and then initialise seq table
  MPI_Bcast(metadata,2,MPI_INT,0,MPI_COMM_WORLD);

  *num_seqs=metadata[0];
  data_size=metadata[1];
  
  //printf("Slave %d knows there are %d sequences of %d words in total\n",
  //	 myid,*num_seqs,data_size);
  seq = (SeqPtr *)  calloc(*num_seqs,sizeof(SeqPtr));
  if (seq == NULL) {
    perror("HandleMPISlaveSetUp allocating memory for seq");
    exit(1);
  }
  seqInfo = (SeqInfoPtr)  calloc(*num_seqs,sizeof(SeqInfoStruct));
  if (seqInfo == NULL) {
    perror("HandleMPISlaveSetUp allocating memory for seqInfo");
    exit(1);
  }
  tree = (UnionFindPtr)  malloc(*num_seqs*sizeof(UnionFindStruct));
  if (tree == NULL) {
    perror("HandleMPISlaveSetUp allocating memory for tree");
    exit(1);
  }
  data= (SeqPtr)  calloc(data_size,sizeof(SeqElt));
  if (data == NULL) {
    perror("HandleMPISlaveSetUp allocating memory for data");
    exit(1);
  }

  maxlen = MPIBUFRECSZ*(*num_seqs);
  buffer = (int32_t *) calloc(maxlen,sizeof(int32_t));
  if (buffer == NULL) {
    perror("Slave buffer allocating memory");
    exit(1);
  }

  // Receive seqInfo stuff
  mpierr(MPI_Bcast(seqInfo,*num_seqs*2,MPI_INT,0,MPI_COMM_WORLD));
  //printf("Slave %d has received sequence metadata\n",myid);

  mpierr(MPI_Bcast(data,data_size,MPI_SHORT,0,MPI_COMM_WORLD));
  //printf("Slave %d has received sequence data <%d>\n",myid,data_size*sizeof(SeqElt));
  recalibrateSeqStruct();
  //printf("Slave %d has recalibrated\n",myid);
  
#endif
}


/* Master distributes work to slaves */
void do_MPImaster_cluster(WorkPtr work) {
#ifdef MPI
  FILE * checkfile;
  int i,k,nlen,tranche,client,maxlen,err,rlen,maxload,minload;
  int bound[2];
  int checkpoint;
  int round=0;
  int  w=(SEQELTWIDTH/2);
  int  *last_sent, *last_got;
  MPI_Status status;
  int seq_data[2];

  maxlen = MPIBUFRECSZ*num_seqs;
  buffer =  (int32_t *) calloc(maxlen,sizeof(int32_t));  
  last_sent = (int *) calloc(numprocs,sizeof(int));
  last_got  = (int *) calloc(numprocs,sizeof(int));
  bzero(last_sent,sizeof(int)*numprocs);
  bzero(last_got,sizeof(int)*numprocs);
  seq_data[0]=num_seqs;
  seq_data[1]=data_size;

  mpierr(MPI_Bcast(seq_data,2,MPI_INT,0,MPI_COMM_WORLD));

  mpierr(MPI_Bcast(seqInfo,2*num_seqs,MPI_INT,0,MPI_COMM_WORLD));
  mpierr(MPI_Bcast(data,data_size,MPI_SHORT,0,MPI_COMM_WORLD));


  /* divide work up */
  tranche = (num_seqs+8)/16/numprocs;
  /* Now wait for requests from slaves */
  for(i=0; i<num_seqs; i=i+tranche) {
    round++;
    if (round < prog_opts.restore) continue;
    checkpoint = 0;
    bound[0]=i;
    bound[1]=MIN(num_seqs,i+tranche);
    if (prog_opts.checkpoint) bound[1]=-bound[1];

    //printf("Master waits for client answer <%d,%d>\n",bound[0],bound[1]);
    merr =MPI_Recv(&client, 1, 
		   MPI_INT, MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD,&status); 
    mpierr(merr);
    //printf("Master gets note from client %d\n",client);
    if (client < 0) { // client wants to send checkpoint 
      client = -client;
      checkpoint = 1;
    }
    // The round previously sent to the client is done
    last_got[client] = last_sent[client];
    // Record the current round send to the client
    last_sent[client]=round;
    // DBg printf("Master sends new work to client %d\n",client);
    merr =MPI_Send(bound, 2, MPI_INT, client, WORKTAG, MPI_COMM_WORLD);
    mpierr(merr);
    if (prog_opts.checkpoint) {
      if (checkpoint) { 
	printf("Master to receive checkpoint\n");
	err = MPI_Recv(buffer, maxlen, MPI_INT, MPI_ANY_SOURCE, 
		       ANSTAG, MPI_COMM_WORLD, &status);
	mpierr(err);
	MPI_Get_count(&status, MPI_INT, &rlen);
	MergeSlaveClusterTable(buffer,rlen);
	checkfile = fopen(prog_opts.checkpoint,"w");
        for(k=1; k<numprocs; k++)
	  fprintf(checkfile,"%d\n",last_got[k]);
	show_clusters(checkfile);
	fclose(checkfile);
      }
      printf("Slave %d sent tranche %d of %d\n",
	     client,round,(num_seqs+1)/tranche);
    }
  }
  bound[0]=-1;
  for(i=1; i< numprocs; i++) {
    // Tell them no more work
    merr = 
      MPI_Recv(&client, 1, MPI_INT, MPI_ANY_SOURCE, 
	       WORKTAG, MPI_COMM_WORLD,&status); 
    //printf("Client %d told to finish\n",client);
    mpierr(merr);
    if (client < 0)  
      client = -client;   
    mpierr(MPI_Send(bound, 2, MPI_INT, client, WORKTAG, MPI_COMM_WORLD));
    err = MPI_Recv(buffer, maxlen, MPI_INT, MPI_ANY_SOURCE, 
		   ANSTAG, MPI_COMM_WORLD, &status);
    mpierr(err);
    MPI_Get_count(&status, MPI_INT, &rlen);
    MergeSlaveClusterTable(buffer,rlen);
  } 
#endif
}








// get work load for the process
// If MPI is not being used then this routine gets the
// workload from the work struct directly. If MPI is
// being used then we send a message to the MPI master
// and return the data
// Where MPI is being used mpi_get_task will be called repeatedly.
// where not, it will only be called once.
int mpi_get_task() {
  int err,bounds[2];

  static int morework=1;
  static int round=0;
  int    fakeid;

  fakeid = -myid;

  if (!morework) return 0;

  if (numprocs==1) { 
    morework=0;  
    return 1;
  }

#ifdef MPI
  MPI_Status status;
  //printf("Client %d Round %d\n",myid,round);
  if ((round%5==2)) 
    MPI_Send(&fakeid, 1, MPI_INT, MASTER, WORKTAG, MPI_COMM_WORLD); 
  else
    MPI_Send(&myid, 1, MPI_INT, MASTER, WORKTAG, MPI_COMM_WORLD); 
  err = MPI_Recv(bounds,2, MPI_INT, MPI_ANY_SOURCE,
		 WORKTAG, MPI_COMM_WORLD, &status);
  //printf("Slave %d got bounds <%d,%d>\n",myid,bounds[0],bounds[1]);
  mpierr(err);
  if (bounds[0]<0) {
    morework=0;
    return 0;
  }
  if (bounds[1]<0) {
    bounds[1]=-bounds[1];
    if (round % 5==2) transmitMPISlaveResponse();
  }

  global_i_beg=bounds[0];
  global_i_end=bounds[1];
  global_j_beg=global_i_beg;
  global_j_end=num_seqs;
  round++;
#endif
  return 1;

}



void wcd_mpi_initialise(int argc, char *argv[]) {
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
}


void wcd_mpi_cleanup() {
#ifdef MPI

      MPI_Finalize();
#endif
}

