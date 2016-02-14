
#ifndef __COMMONINCL
#define __COMMONINCL

#include <stdio.h>

/* For FreeBSD, you may have to replace <stdint.h> with  <sys/types.h> */
#include <stdint.h>


#include <stdlib.h>

#define DEBUGASSERT 
#ifdef DEBUGASSERT
#define ASSERT(x) assert(x);
#else
#define ASSERT(x)
#endif


/* The following are constants that the user might want to change */

/* Longest EST we'll see is */
#define MAX_SEQ_LEN 512000
#define MAX(i,j) (i>j?i:j)
#define MIN(i,j) (i<j?i:j)


#ifdef WIN32
#define HIST_MAX 100000
#else
#define HIST_MAX 1000000
#endif

/* used for num_common_matches */

/* BIGWORDTSIZE == 4^sample_word_len */
#define SAMPLE_WORD_LEN 8
#define BIG_WORD_TSIZE (1<<(SAMPLE_WORD_LEN<<1))
#define BIG_WORD_MASK ((BIG_WORD_TSIZE-1))
#define BIG_TOP_MASK (3*(BIG_WORD_TSIZE>>2))
#define BIG_SHIFT ((SAMPLE_WORD_LEN-1)<<1)
#define GAP 6


#define VBIG_WORD_LEN 13
// BIGWORDTSIZE == 4^VBIG_WORD_LEN
#define VBIG_WORD_TSIZE (1<<(VBIG_WORD_LEN<<1))
#define VBIG_WORD_MASK ((VBIG_WORD_TSIZE-1))
#define VBIG_TOP_MASK (3*(VBIG_WORD_TSIZE>>2))
#define VBIG_SHIFT ((VBIG_WORD_LEN-1)<<1)

#define IGNORE  1
#define FIX     2
#define RESET   5
#define DUSTED  7
//nb: reset implies not ignore

#define IS_FLAG(flagvalue,x)   ((seqInfo[x].flag&flagvalue) == flagvalue)

#define SET_FLAG(flagvalue,x)  (seqInfo[x].flag ^= flagvalue)
#define IGNORE_SEQ(x)           ((IS_FLAG(IGNORE,x)) != ignore_active)

#define CODE(arg)  (tr[((arg>>1)&3)])

//NB: sizeof(wordElt) > sample_word_len

typedef uint16_t wordElt; 
typedef unsigned short SeqElt;


typedef SeqElt  * SeqPtr;

#define  SEQELTWIDTH  (sizeof(SeqElt)*8)
#define  BASESPERELT  (SEQELTWIDTH/2)
#define  SEQBYTESIZE  (SEQELTWIDTH/8)
 


typedef char * string;


typedef struct UnionFindStruct {
  int     cluster; // refers to index of root
  int     next;    // for find-union: next in linked list of cluster
  int     last;    // for find-union: last ditto
  int     match;   // index of another seq for which a match was made
  short   orient;  // orientation of that sequence
  short   rank;    // for find-union   
} UnionFindStruct;


typedef UnionFindStruct *UnionFindPtr;

typedef struct SeqInfoStruct  {
  // NB len and flag *must* be next to each other
  // -- MPI transmits the two together
  int     len;     // length
  int     flag;    // flags about status
} SeqInfoStruct;



typedef SeqInfoStruct * SeqInfoPtr;


typedef struct SeqIDStruct {
  string  id;      // sequence id (from data file)
  string  clone;   // clone id
} SeqIDStruct;


/* NB: I split the UF and SeqStruct. The former is dynamic -- changes
   during computation, whereas latter doesn't change after the data is
   read in. I experimented with having each WorkBlock have its own UF
   structure. This would eliminate the need for synchronisation at the
   cost of some extra memory and passing lots more parameters around.
   Experimentation showed no difference. The reason for poor pthreads
   performance is therefore not due to synchronisation costs. I therefore
   restored only having one UF structure, but it makes sense to keep the
   UF and SeqStructs separate for future changes. Also, there's a small 
   MPI savings to be made */

typedef SeqIDStruct   * SeqIDPtr;



//----------------- The sequences



// w contains a word. Now we want to shift another char into the word
// and shift the oldest char out
//   seq -- sequence from where we get the char
//     i -- which char it is in the sequence
//           NB: we store the sequence in compact form so 
///          it's a little more complex than seq[i]

//#define GETWORD(w,seq,i,shf,mask) ((((seq[i/BASESPERELT]>>((i%BASESPERELT)*2))&3)|(w<<2))&mask)

#define GETWORDRC(w,seq,i,shf,mask)   (((((~seq[i/BASESPERELT])>>((i%BASESPERELT)*2))&3)|(w<<2))&mask)

#define GETWORD(w,seq,i,shf,mask) (((((seq[i/BASESPERELT]>>((i%BASESPERELT)*2))&3)<<shf)|(w>>2))&mask)






#ifdef PTHREADS
#define PLOCK(y)    pthread_mutex_lock(y); 
#define PUNLOCK(y)  pthread_mutex_unlock(y); 
#else
#define PLOCK(y) 
#define PUNLOCK(y) 
#endif
          


/* Only change what's below here if you know what you're doing
   -- warning doesn't apply to original author   */


// These are 'global' variables. Some could be local but are
// declared global to avoid repeated allocation and deallocation


typedef struct WorkBlock {
  int thread_num;
  int i_beg;
  int i_end;
  int j_beg;
  int j_end;
  short *tableP;
  short *tableR;
  char    tableQ[VBIG_WORD_TSIZE];
  int     word_seen[MAX_SEQ_LEN]; 
  wordElt word_seenQ[MAX_SEQ_LEN];
  wordElt * word1;
  wordElt * word2;
  void * fn_data;
  char * filename;
  int    do_dump;
  FILE * dump_file;
  int    workflag;
  int  * index;
} WorkBlock;

typedef WorkBlock * WorkPtr;

typedef int (*distFunType)(WorkPtr work, int s1, int s2, int rc);

distFunType  dist, distpair;


void fasta(FILE *outf, unsigned int i);

void get_bounds(WorkPtr work, int s1, int s2, int * left, int * right);

void create_word_lists(WorkPtr work, int s1, int s2, int rcflag );

void set_up_word_table(WorkPtr work, int s1);

void clear_word_table(WorkPtr work, int s1);

int find_parent(int i);

int mini_find_parent(int i);

void invert_orientation(int k);

void make_uniond2(int i, int j, int invert);

void show_clusters(FILE * outf);

void show_EXT(FILE *outf);

void num_matches(WorkPtr work, int s1, int s2, int *posmat, int *rcmat);

void stat_num_matches(WorkPtr work, int s1, int s2, int *posmat, int *rcmat);

void read_sequences(FILE *finp, int c1, int c2);


void produce_clusters(int d, char *dirname);

int under_threshold(int num_mat, int threshold);

int similar(WorkPtr work, int i, int j, distFunType dist, int theta);

void initialise(WorkPtr work, FILE * edfile);

char cod2ch(int i);

void init_sequences(FILE *finp, int c1, int c2);

void sample_heuristic(WorkPtr work, int s1, int s2, int *posmat, int
*rcmat);


void tv_heuristic(WorkPtr work, int s1, int s2, int *posmat, int *rcmat);

int tv_heuristic_rc(WorkPtr work, int s1, int s2);

int tv_heuristic_pos(WorkPtr work, int s1, int s2) ;

void debug_show_clusters(WorkPtr work, char * fname);

int ind(int i);

void  mirror_sequences(char *);

int compare_int(const void *a, const void *b);

void pseq(int i);
void pseqi(int i);


#endif
