
/* 	$Id$	 */


#define WORDS_BIGENDIAN

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
#include <sys/mman.h>
#import "math.h"


#define SA_LEN 10
#define lineBlockSize 32
#define TABLESIZE 512000
#define HIST_MAX  100000
int histo [HIST_MAX];
int maxhist=0;

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

int size1, * line;
char * sequences;
int *sa, *inv, * cout, * lftind, * rgtind, *count;
int seqFileSize,sarrayFileSize;
int *candidates, *matches;
int indices;
int glob_bmax=0;




#define BUFFLEN 8000
#define BUFFSIZE   (BUFFLEN*4)


int read_int() {  
  // read the next integer from the suffix array file
  char conv[4];
  int ans;
  static char buff[BUFFSIZE]; // keep a buffer
  static int k=BUFFSIZE;             // index into buffer
  static int n=BUFFSIZE;       // num bytes read in from buffer
  // check  whether we have used buffer or if empty 
  if (k>=n) { // buffer must be exhausted
    // check if last time we read in we hit the end of the file
    if (n<BUFFSIZE) return -1; 
    // now read in the buffer
    n=read(indices,buff,BUFFSIZE);
    if (n<4) return -1;
    k=0;
  }
#ifdef WORDS_BIGENDIAN
  conv[3]=buff[k+3];
  conv[2]=buff[k+2];
  conv[1]=buff[k+1];
  conv[0]=buff[k+0];
#else
  conv[0]=buff[k+3];
  conv[1]=buff[k+2];
  conv[2]=buff[k+1];
  conv[3]=buff[k+0];
#endif
  k=k+4;
  ans = *((int *) conv);
  return ans;
}





int get_indices(char * sarray_name,int * sa, int * inv, int * line, int size) {
  // open suffix array and read it and create inverse
  // returns biggest # times a word is repeated
  int skip,i,k, seqnum, curr_word, curr_ind,bsize,bmax;

  char  oldstr[64];
  // These vars are used only in stats analysis
#ifdef WORDSTATS
  float eta;
  int   num_unique_words=0,bsize,total=0;
  int   table[32];
  memset(table,0,32*sizeof(int));
#endif  
 
  // read in suffix array
  indices = open(sarray_name,0,0);
  if (indices<0) {
    perror(sarray_name);
    exit(1);
  }
    
  for(i=0; i<size; i++) {
    k = read_int();
    sa[i]=k;
    inv[k]=i;
  }  
  close(indices);


  // we can skip first part
  sprintf(oldstr,"%d",num_seqs);
  skip = 0;//2*num_seqs*(strlen(oldstr)+1)-20;

  curr_word = sa[skip];
  i=curr_ind  = skip;
  while(i<size) {
    while((i < size) && (strncmp((char *) (sequences+curr_word),(char *) (sequences+sa[i]),suffix_len)==0)) {
      //assert(sa[i]<size);
      inv[sa[i]]=curr_ind;
      i++;
    }
    bsize = i-curr_ind;
    if (bsize>bmax) bmax=bsize;
#ifdef WORDSTATS
    total = total + bsize;
    num_unique_words++;
#endif
    curr_ind = i;
    curr_word = sa[i];

  }
  

#ifdef WORDSTATS  
  eta=total*1.0/num_unique_words;
  printf("%d %d %d %f\n",suffix_len,num_unique_words,bmax,eta);

  for(i=0; i<16; i++)
    printf("%5d-%5d: %10d\n",(1<<i),(2<<i)-1,table[i]);
#endif

  seqnum=0;
  for(i=0; i<size; i++) {
    line[i]=seqnum%(num_seqs/2);
    if ((int)  sequences[i]==-1) seqnum++;
  }
  
  return bmax;
 
}


int scan_next_seq(int i, int q, int * num_mat)  {
  int start_match, curr_match,num_block, j,z=0;
	
  while ((int) sequences[q] != -1) {
    start_match = curr_match = inv[q];  // curr_match points to sa
    // j  is a sequence id
    //assert(inv[sa[curr_match]]==curr_match);
    //printf("Sequence %d : %s\n",curr_seq,word);
		
    num_block=0;
    while (inv[sa[curr_match]]==start_match) {
      j = line[sa[curr_match]];
      curr_match++;
      if (j<i) continue;
      if (lftind[j]<0) {
	matches[z]=j;
	z++;
	lftind[j]=rgtind[j]=q;
      } else
	rgtind[j]=q;
      count[j]++;
      num_block++;
    }
    q++;
  }
  * num_mat=z;
  return q+1;
}


inline
int check_candidates(int i, int num_mat) {
  int r,j;
  // now check for each matching sequence see if good enough match
  int num_cand=0;
  for(r=0; r<num_mat; r++) {
    j = matches[r];
    if ((j>i) && (count[j]>=alpha) && (rgtind[j]-lftind[j]>=beta)) {
      candidates[num_cand]=j%(num_seqs/2);
      num_cand++;
    }
    count[j] = 0;
    lftind[j]=-1;
  }
  return num_cand;
	
}

inline
void scan_seqs_for_words(WorkPtr work,int * sa, int * inv,  int * line, int size) {

  int curr_seq=0, i=0, q=0, num_mat=0, j;
  int num_cand=0;
  char word[64];
	
  memset(word,0,60);

  while(q<size/2-suffix_len) {
    i = line[q];
    //if (i%100==0) printf("Seq %d\n",i);
    //assert(i == curr_seq);          // q points to data

    // work through all words in current sequence
    q=scan_next_seq(i,q,&num_mat);
    // find all very good matches
    num_cand = check_candidates(i,num_mat);

    /*
    histo[num_cand]++;
    if (num_cand>maxhist) maxhist=num_cand;
    if (num_cand>3000) {
      printf("%d",i);
      for(j=0; j<num_cand; j++)
         printf(" %d",candidates[j]);
      printf("\n");
    }
    */
    // for each good match cluster
    if (num_cand>=1) {
      //printf("Cluster %d: candidates %d\n",i,num_cand);
      complete_pairwise_cluster(work,i,candidates,num_cand);
    }
		
    // Now skip the RC -- it will have been picked up already

    //q = q + (int) (0.9*seqInfo[i].len);
    //while(line[q]>i) q--;
    //q--;
    //while ((int) sequences[q] != -1) q++;
    curr_seq++;
    //q++;

  }

}
			
		

void do_suffix_cluster (WorkPtr work) {

  struct stat st;
  int fd;
  int size0, bmax;
  
  char nlcname[128], sarray_name[128];

  printf("There are %d seqs\n",num_seqs);
  lftind = (int *) calloc(num_seqs,sizeof(int));
  rgtind = (int *) calloc(num_seqs,sizeof(int));
  count  = (int *) calloc(num_seqs,sizeof(int));

  memset(histo,0,HIST_MAX*sizeof(int));
  memset(count,0,num_seqs*sizeof(int));
  memset(lftind,-1,num_seqs*sizeof(int));
  sprintf(nlcname,"%s.%s",work->filename,"ois");
  sprintf(sarray_name,"%s.%s",work->filename,"suf");

  stat(nlcname, &st);
  seqFileSize = st.st_size;
  fd = open(nlcname, 0, 0);
  if (fd == -1) {
    perror(nlcname);
    exit(2);
  }
  size0 = st.st_size;
  size1 = getpagesize() * ((size0 + getpagesize()) / getpagesize());
  sequences =
    mmap((void *) 0, size1, PROT_READ | PROT_WRITE,
	 MAP_PRIVATE, fd, 0);
  if (sequences == MAP_FAILED) {
    printf("mmap() failed on the data file\n");
    exit(2);
  }
  sequences[size0]=-1;


  sa   = (int *) calloc(size0, sizeof(int));   // elts at least log(mn) bits
  inv  = (int *) calloc(size0, sizeof(int));
  line = (int*)  calloc(size0, sizeof(int));  // elts at least log(n) bits
 
  bmax = get_indices(sarray_name,sa,inv,line,size0);
  
  candidates = (int *) calloc(num_seqs,sizeof(int));
  matches = (int *) calloc(num_seqs+16,sizeof(int));

  printf("==================================\n");
  scan_seqs_for_words(work,sa,inv,line,size0);

  num_seqs=num_seqs/2;


}

void  suffix(WorkPtr work, int s1, int s2, int rcflag) {
  // dummy -- this is nevery called but we do use it
  // in wcd.c to assign to do_cluster. We can then test
  // do_cluster to check what to do
  return;
}
