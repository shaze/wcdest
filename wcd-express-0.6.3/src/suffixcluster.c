    
/* 	$Id$	 */

typedef long long unsigned s_int;

#define SUFFSIZE (sizeof(s_int))
#define MAXSUFFBIT (SUFFSIZE-1)

int SKIPPED=0;

#undef SYMMETRIC
#undef WORDSTATS


#ifndef lint
static char vcid[] = "$Id$";
#endif /* lint */


#include "common.h"
#include "wcd.h"
#include "d2.h"

#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#ifndef WIN32
#include <sys/mman.h>
#else
#include <windows.h>
#endif
#import "math.h"


#define BLOCKFLAG (1<<30)
#define BLOCKMASK (BLOCKFLAG-1)
#define SA_LEN 10
#define lineBlockSize 32

//#define HIST_MAX  100000
//
//int histo [HIST_MAX];
//int maxhist=0;

extern int boost;
extern int ignore_active;
extern int word_len;
extern int window_len;
extern int num_seqs;
extern int alpha;
extern int beta;

extern SeqPtr       data;
extern SeqInfoPtr   seqInfo;
extern UnionFindPtr tree;
extern int theta;
extern int suffix_len;
extern int NUM_num_matches, total_matches, NUM_dfn;

s_int size0, size1;
int   * line;        // for each byte in data, which seq it belongs to
char  * sequences;   // pointer to all data
s_int * seqtextptr;  // for each sequence, index to where in sequences 
s_int *sa, *inv;     // suffix array and its invoice

int  * lftind, * rgtind;
int *count;
#ifdef SYMMETRIC
int *j_lftind, *j_rgtind;
#endif

int sarrayFileSize;
int *candidates, *matches;
int indices;
int glob_bmax=0;


#ifdef WORDSTATS
unsigned long long int eta=0,etasq=0;
long int   num_unique_words=0,bsize,total=0;
int   totcands=0;
int   table[32];
#endif  



#define BUFFLEN 8000
#define BUFFSIZE   (BUFFLEN*4)


s_int read_int() {  
  // read the next integer from the suffix array file
  char        conv[SUFFSIZE];
  s_int       ans;
  int         i;
  static char buff[BUFFSIZE]; // keep a buffer
  static int  k=BUFFSIZE;             // index into buffer
  static int  n=BUFFSIZE;       // num bytes read in from buffer

  // check  whether we have used buffer or if empty 
  if (k>=n) { // buffer must be exhausted
    // check if last time we read in we hit the end of the file
    if (n<BUFFSIZE) return -1; 
    // now read in the buffer
    n=read(indices,buff,BUFFSIZE);
    if (n<SUFFSIZE) return -1;
    k=0;
  }
#ifdef WORDS_BIGENDIAN
  for(i=0; i<SUFFSIZE; i++)
    conv[i]=buff[k+MAXSUFFBIT-i];
#else
  for(i=0; i<SUFFSIZE; i++)
    conv[i]=buff[k+i];
#endif
  k=k+SUFFSIZE;
  ans = *((s_int *) conv);
  return ans;
}



int compare_addr(const void *a, const void *b) {
  int diff = line[*((s_int *) a)]-line[*((s_int *) b)];
  return diff;
}



void sort_suffix_array(s_int * sa, s_int size, int flag) {
  // sort suffix array on sequence position
  // also creates the inverse of the suffix array
  // -- don't do sort if we are creating symmetric adjacency list
  //    as output
  s_int curr_word,curr_ind,i,bsize,k;
  char  www[255];

  curr_word = sa[0];
  i=curr_ind=0;
  // Find each block and then sort each block separately
  while(i<size) {
    while((i < size) && 
          (strncmp((char *) (sequences+curr_word),
		   (char *) (sequences+sa[i]),suffix_len)==0)) {
      i++;
    }

    bsize = i-curr_ind;
#ifdef WORDSTATS
    if (bsize>glob_bmax) glob_bmax=bsize;
    eta   = eta+bsize;
    etasq = etasq+((long long) bsize)*bsize;
    num_unique_words++;
#endif
    if ((bsize>1)&&(flag!=GO_FLAG_ADJSYMM))
      qsort(&sa[curr_ind],(int) bsize,sizeof(s_int),compare_addr);

    line[sa[curr_ind]] = line[sa[curr_ind]] | BLOCKFLAG;
    curr_ind = i;
    curr_word = sa[i];
  }
  // This handles the last block
  bsize = i-curr_ind;
  if (flag != GO_FLAG_ADJSYMM) 
    qsort(&sa[curr_ind],(int)bsize,sizeof(s_int),compare_addr);
  for(i=0; i<size; i++) {
      if (sa[i]<size/2)
       inv[sa[i]]=i;
  }

  line[sa[i]] = line[sa[i]] | BLOCKFLAG;
}


int get_indices(char * sarray_name,s_int * sa, s_int * inv, int * line, 
                s_int size, int flag) {
  // open suffix array and read it and create inverse
  // returns biggest # times a word is repeated
  int j, seqnum, w1_s, w2_s, w1_f, w2_f;
  s_int i,k;

  seqnum=0;
  seqtextptr[0]=0;
  for(i=0; i<size; i++) {
    line[i]=seqnum%(num_seqs/2);
    if ((int)  sequences[i]==-1) seqnum++;
  }

  // read in suffix array
  indices = open(sarray_name,0,0);
  if (indices<0) {
    perror(sarray_name);
    exit(1);
  }
    
  for(i=0; i<size; i++) {
    k = read_int();
    sa[i]=k;
  }  
  close(indices);

  sort_suffix_array(sa,size,flag);

  k=j=0;

  seqnum=0;

  /*
  for(i=0; i<suffix_len; i++) 
    if (sequences[i]=='X') k++;

  */
  for(i=suffix_len; i<size/2-suffix_len; i++) {
    if ((int)  sequences[i]==-1) {
      j=0;
      seqnum++;
      seqtextptr[seqnum]=i+1;
    }

    //if (sequences[i]=='X')  {j++; k++;}
    //if (sequences[i-suffix_len]=='X') k--;
    //if (k>=suffix_len/3) inv[i-suffix_len] = size1;

  }

  /*
  j=k=0;
  for(i=0; i<window_len; i++) 
    if (sequences[i]=='X') k++;
  for(i=window_len-suffix_len; i<2*window_len-suffix_len; i++)
    if (sequences[i]=='X') j++;
  w1_s=0;
  w1_f=window_len-1;
  w2_s=window_len-suffix_len;
  w2_f=w2_s+window_len-1;
  for(i=window_len-suffix_len; i<size1/2-suffix_len; i++) {
    if ((j>window_len/4) && (k>window_len/4))
      inv[i]=size1;
    if (sequences[w1_s]=='X') k--;
    w1_s++;
    w1_f++;
    if (sequences[w1_f]=='X') k++;
    if (sequences[w2_s]=='X') j--;
    w2_s++;
    w2_f++;
    if (sequences[w2_f]=='X') j++;
  }
  */
  
  return 0;
 
  }




s_int scan_next_seq(int i, int * num_mat)  {
  int  j,z=0,bsize;
  s_int where, start_match, curr_match, q, start;	

  start = seqtextptr[i];
  if (IGNORE_SEQ(i)) {
    return start+seqInfo[i].len-1;
  }
  for(q=start; sequences[q] != -1; q++) {
    start_match = curr_match = inv[q];  // curr_match points to sa
    // j  is a sequence id
    /*
    if (start_match==size1) {
      SKIPPED++;
      continue;
    }
    */
    while (1) {
      where=sa[curr_match];
      j = line[where] & BLOCKMASK;
      /* Could check if j should be ignored too -- but this will
       * be relatively infrequently -- check in candidates rather */
      curr_match++;
      if (lftind[j]<0) {
	matches[z]=j;
	z++;
	lftind[j]=rgtind[j]=q-start;
#ifdef SYMMETRIC
	j_lftind[j]=j_rgtind[j]=where-starttextptr[j];
#endif	
      } else {
	rgtind[j]=q-start;
#ifdef SYMMETRIC
	j_lftind[j]=MIN(j_lftind[j],where-starttextptr[j]);
	j_rgtind[j]=MAX(j_rgtind[j],where-starttextptr[j]);
#endif
      }
      count[j]++;
      if (line[sa[curr_match]] & BLOCKFLAG) break;
    }
    bsize = curr_match-start_match+1;
  }
  *num_mat=z;
  return q;
}



int kl_check_candidates(int i, int num_mat) {
  int r,j;
  // now check for each matching sequence see if good enough match
  int num_cand=0;
  for(r=0; r<num_mat; r++) {
    j = matches[r];
    if ((count[j]>=alpha) && (rgtind[j]-lftind[j]>=beta)
#ifdef SYMMETRIC
        && (j_rgtind[j]-j_lftind[j]>=beta)
#endif
    ) {
      ASSERT(j<num_seqs);
      candidates[num_cand]=j%(num_seqs/2);
      num_cand++;
    }
    count[j] = 0;
    lftind[j]=-1;
  }
  ASSERT(num_cand<num_seqs);

  return num_cand;
	
}



int check_candidates(int i, int num_mat) {
  int r,j,a,b;
  // now check for each matching sequence see if good enough match
  int num_cand=0, pos;
  if (seqInfo[i].len<window_len) {
    a=1;
    b=beta>>2;
  } else {
    a=alpha;
    b=beta;
  }
  for(r=0; r<num_mat; r++) {
    j = matches[r];
    if (IGNORE_SEQ(j%(num_seqs/2))) continue;
    pos = count[j]>=a;
    if ((j>i) && (pos) && (rgtind[j]-lftind[j]>=b)
#ifdef SYMMETRIC
        && (j_rgtind[j]-j_lftind[j]>=b)
#endif
    ) {
      ASSERT(j<num_seqs);
      candidates[num_cand]=j%(num_seqs/2);
      num_cand++;
    }
    count[j] = 0;
    lftind[j]=-1;
  }
#ifdef WORDSTATS
  totcands = totcands+num_cand;
#endif
  ASSERT(num_cand<num_seqs);

  return num_cand;
	
}


void scan_seqs_for_links(FILE * outf, WorkPtr work,
			 s_int * sa, s_int * inv,  int * line, s_int size) {

  int curr_seq=0, i=0, num_mat=0, j;
  s_int q;
  int num_cand=0;

  for(i=0; i<num_seqs/2; i++) {
    // work through all words in current sequence
    q=scan_next_seq(i,&num_mat);
    // find all very good matches
    num_cand = kl_check_candidates(i,num_mat);

    // for each good match cluster
    if (num_cand>=2) {
      fprintf(outf,"%d:",i);
      qsort(candidates,num_cand,sizeof(int),compare_int);
      complete_klink_prep(outf,work,i,candidates,num_cand);
    }
    else
      fprintf(outf,"%d:%d.\n",i,i);
    curr_seq++;

  }

}


void scan_seqs_for_words(FILE * outf, WorkPtr work,
			 s_int * sa, s_int * inv,  
			 int * line, int size) {

  int curr_seq=0, i=0, num_mat=0, j;
  s_int q=0;
  int num_cand=0;

  for(i=0; i<num_seqs/2; i++) {
    q= scan_next_seq(i,&num_mat);
    // find all very good matches
    num_cand = check_candidates(i,num_mat);

    // for each good match cluster
    if (num_cand >=1) 
      complete_pairwise_cluster(work,i,candidates,num_cand);

    curr_seq++;
  }

}
			

void i_do_suffix_cluster (FILE * outf, WorkPtr work) {

  struct stat st;
  int fd, bmax;
  
  char nlcname[512], sarray_name[512];

  // Allocate memory and initialise
  lftind = (int *) calloc(num_seqs,sizeof(int));
  rgtind = (int *) calloc(num_seqs,sizeof(int));
  count  = (int *) calloc(num_seqs,sizeof(int));
  seqtextptr = (s_int *) calloc(num_seqs,sizeof(s_int));
#ifdef SYMMETRIC
  j_lftind = (int *) calloc(num_seqs,sizeof(int));
  j_rgtind = (int *) calloc(num_seqs,sizeof(int));
#endif

#ifdef WORDSTATS
  memset(table,0,32*sizeof(int));
#endif

  memset(count,0,num_seqs*sizeof(int));
  for (fd=0; fd<num_seqs; fd++) 
    lftind[fd]=-1;

  // Open suffix array  and read in
  sprintf(nlcname,"%s.%s",work->filename,"ois");     //Stripped data
  sprintf(sarray_name,"%s.%s",work->filename,"suf"); //Suffixes

  stat(nlcname, &st);
  size0 = st.st_size;
  fd = open(nlcname, 0, 0);
  if (fd == -1) {
    perror(nlcname);
    exit(2);
  }
  size1 = getpagesize() * ((size0 + getpagesize()) / getpagesize());
#ifndef WIN32
  sequences =
    mmap((void *) 0, size1, PROT_READ | PROT_WRITE,
	 MAP_PRIVATE, fd, 0);
  if (sequences == MAP_FAILED) {
    printf("mmap() failed on the data file\n");
    exit(2);
  }
#endif
  sequences[size0]=-1;


  sa   = (s_int *) calloc(size0+sizeof(s_int), sizeof(s_int));   
  inv  = (s_int *) calloc(size1/2, sizeof(s_int));
  line = (int*)    calloc(size0+sizeof(s_int), sizeof(int));  
 
  get_indices(sarray_name,sa,inv,line,size0,work->workflag);
  
  candidates = (int *) calloc(num_seqs,sizeof(int));
  matches    = (int *) calloc(num_seqs+16,sizeof(int));

  if (outf == NULL)
    scan_seqs_for_words(outf,work,sa,inv,line,size0);
  else
    scan_seqs_for_links(outf,work,sa,inv,line,size0);

  num_seqs=num_seqs/2;

#ifdef WORDSTATS
   printf("%s suffixlen %d numseqs %d uniquewords %d eta %ld %6.4f  etasq %ld %6.4f etamax %d totcands %d NUMM %d D2calls %d d2succ %d\n ",
	  work->filename,
	  suffix_len,num_seqs,num_unique_words,
	  eta,1.0*eta/num_unique_words, etasq,1.0*etasq/num_unique_words,glob_bmax,totcands,
      NUM_num_matches,NUM_dfn,total_matches);
  printf("etaseq=%ld, gamma=%5.3f etasq=%6.3f\n",etasq,1.0*num_unique_words/size0,1.0*etasq/num_unique_words);
   printf("SKIPPED %d\n",SKIPPED);
#endif


}



void do_suffix_cluster (WorkPtr work) {
  FILE *outf=NULL;
  i_do_suffix_cluster(outf,work);
}

void do_kseed_suffixcluster(FILE *outf, WorkPtr work) {
  i_do_suffix_cluster(outf,work);
}


void  suffix(WorkPtr work, int s1, int s2, int rcflag) {
  // dummy -- this is nevery called but we do use it
  // in wcd.c to assign to do_cluster. We can then test
  // do_cluster to check what to do
  return;
}
