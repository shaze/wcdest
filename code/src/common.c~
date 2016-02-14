


#ifndef lint

     #endif /* lint */



#ifndef __WCDHINCL
#include "wcd.h"
#endif

#ifndef __D2HINCL
#include "d2.h"
#endif

#ifndef __EDHINCL
#include "ed.h"
#endif

#ifndef __COMMONINCL
#include "common.h"
#endif


#include "suffixcluster.h"

#ifdef PTHREADS
#include "pthread.h"
extern pthread_mutex_t write_union_mutex, find_parent_mutex;
#endif

#include <strings.h>
#include <string.h>


#ifdef WIN32
#define random rand
#define srandom srand
#endif


// -- global vars used 
int  tr[4];  // see note below -- for translating input
SeqElt *rc, *rc_big; // table for RC

//----------  parameters for clustering
int alpha=3;
int beta=45;
int boost=0;
int clonelink = 0;
int num_threads=1;
int num_words_win;      // number of words in window
int rc_check = 1;
int sample_thresh=4;
int sample_thresh_def;
int skip=1;
int suffix_len = 16;
int theta = 40; // default value
int theta_def;
int window_len = 100;
int window_len_def;
int word_len   = 6;
int word_mask;   // word_tszie1
int word_shift;  // 2*(word_len-1)
int word_threshold=70;
int word_threshold_def;
int word_tsize;  // 4^word_len
//---------- flags for output and behaviour etc
int sp, sr;
int ignore_active=0;
int reindex_value=0;

//--- the following are primarily for stats
int NUM_num_matches = 0;
int NUM_dfn = 0;
int NUM_dfn_succ = 0;
int NUM_full_heuristic=0;
int total_matches = 0;
int NUM_pairwise=0;
FILE * outf;


//-------------- This is the main sequence array


SeqPtr  *seq;
SeqPtr  data;

SeqInfoPtr  seqInfo;
UnionFindPtr   tree;
SeqIDPtr    seqID;

int num_seqs=0;  // how many there are


//-----------

int myid = 0;
int numprocs=1;








//----------------------------------------------------
// Given two sequences s1 and s2 create two arrays of the words
// that are in the sequences and say how many words there are
// first some auxiliary functions



int compare_int(const void *a, const void *b) {
  return  (* ((int *) a)) - (* ((int *) b));
}



void set_up_word_table(WorkPtr work, int s1) {
  int i;
  uint16_t w1=0;

  for(i=0; i<SAMPLE_WORD_LEN-1; i++) {
    w1=GETWORD(w1, seq[s1], i,BIG_SHIFT, BIG_WORD_MASK);
  }
  for(i=SAMPLE_WORD_LEN-1; i < seqInfo[s1].len;i++) {
    w1  =  GETWORD(w1,seq[s1],i,BIG_SHIFT, BIG_WORD_MASK);
    work->tableP[w1]=1;
    work->tableR[rc_big[w1]]=1;
    work->word_seen[i]=w1;
  }
}


void clear_word_table(WorkPtr work, int s1) {
  int i;
  for(i=SAMPLE_WORD_LEN-word_len; i<seqInfo[s1].len ; i++) {
    work->tableP[work->word_seen[i]]=0;
    work->tableR[rc_big[work->word_seen[i]]]=0;
  }
}



// union_find algorithms ---------------------------------------




int find_parent(int i) {
  int k, r;

  r = i;

  while (r != tree[r].cluster)
    r = tree[r].cluster;
  k = i;

  i = tree[k].cluster;
  if (i != tree[i].cluster) { 
    while (k != tree[k].cluster) {
      i = k;
      k = tree[k].cluster;
      tree[i].cluster=r;
    }
    ASSERT(tree[tree[r].last].next<0);
  }
  return r;
}




void invert_orientation(int k) {
  // When we take the union of two clusters, we must make sure
  // that the orient field of the cluster that will be the child 
  // of the root is kept in the correct orientation with respect
  // to the root. 
  // If i and j were the sequence matches then
  //   -- if the match was a positive match, the orientation
  //      must change exactly when the orientation of i wrt to
  //      its root is different to the orientation of j wrt to its root.
  //   -- if the match was an rc match, the orientation must change
  //      exactly when the two orientations of i and j are the same.
  int m;
  for (m=k ; m>=0; m=tree[m].next) {
    tree[m].orient = -tree[m].orient;
  }
}


int mini_find_parent(int i) {
  while (i != tree[i].cluster)
    i = tree[i].cluster;
  return i;
}

/* merge two existing unions */

void make_uniond2(int i, int j, int invert) {
  int q;

  PLOCK(&write_union_mutex);
  PLOCK(&find_parent_mutex);
  i = find_parent(i);
  j = find_parent(j);
  if (i == j) {
    PUNLOCK(&find_parent_mutex);
    PUNLOCK(&write_union_mutex);
    return;
  }

  if (tree[i].rank > tree[j].rank) {
#ifdef MPIVERSION
    if (invert) invert_orientation(j);
#endif
    q = tree[i].last;
    tree[q].next = j;
    tree[i].last = tree[j].last;
    tree[j].cluster = i;
    if (IS_FLAG(FIX,j)) SET_FLAG(FIX,i);
    ASSERT(tree[tree[i].last].next<0);
  } else {
#ifdef MPIVERSION
    if (invert) invert_orientation(j);
#endif
    q = tree[j].last;
    tree[q].next = i;
    tree[j].last = tree[i].last;
    tree[i].cluster = j;
    if (tree[i].rank==tree[j].rank) tree[j].rank++;
    if (IS_FLAG(FIX,i)) SET_FLAG(FIX,j);
    ASSERT(tree[tree[j].last].next<0);
  }
    PUNLOCK(&find_parent_mutex);
    PUNLOCK(&write_union_mutex);


}



void show_clusters(FILE *outf) {
  int i, j;
  for (i=reindex_value; i<num_seqs; i++) {
    if (tree[i].cluster == i) {
      fprintf(outf,"%d",i);
      for (j=tree[i].next; j>=0; j = tree[j].next) {
        fprintf(outf," %d",j);
      }
      fprintf(outf,".\n");
    }
  }
}

void debug_show_clusters(WorkPtr work, char * fname) {
  FILE * out;
  out = fopen(fname,"w");
  show_clusters(out);
  fclose(out);
}


void make_rc(SeqElt * rc, int width, int maxn) {
     int i, j, w;
     for (i=0; i<maxn; i++) {
       rc[i]=0;
       w = ~i;
       for(j=0; j<width; j++) {
	 rc[i] = (rc[i] << 2) | (w & 3);
	 w = w>>2;
       }
     }
}
       

uint code(char ch) {
  switch(ch) {
  case   'A' :
  case   'a' : return 0;
  case   'C' :
  case   'c' : return 1;
  case   'G' :
  case   'g' : return 2;
  case   'T' :
  case   't' : return 3;
  default:
    return (random()>>5)&3;
  }
}



char cod2ch(int i) {
  switch (i) {
  case 0 : return 'A'; break;
  case 1 : return 'C'; break;
  case 2 : return 'G'; break;
  case 3 : return 'T'; break;
  default :
    return 'M';
  }
}



void fasta(FILE *outf, uint i) {  
  uint k;
  int  w=(SEQELTWIDTH/2);
#ifndef NOAUXINFO
  fprintf(outf,">%s",seqID[i].id);
#else
  fprintf(outf,">%07d",i);
#endif
  for(k=0; k<seqInfo[i].len; k++) {
    if (k%80==0) fprintf(outf,"\n");
    fprintf(outf,"%c",cod2ch(seq[i][k/w]>>((k%w)*2)&3));
  }
  fprintf(outf,"\n");
}





void produce_clusters(int d, char *dirname) {
  int i, j, k=1, count;
  char fname[255];
  FILE *S0,  *T;
  sprintf(fname,"%s/S0.fasta",dirname);
  S0 = fopen(fname,"w");
  if (S0 == NULL) {
       perror(fname);
       exit(1);
  }
  for (i=reindex_value; i<num_seqs; i++) {
    if (tree[i].cluster == i) {
      count=1;
      for (j=tree[i].next; j>=0; j = tree[j].next) {
        count++;
      }
      if (count<d) 
	T = S0;
      else {
	sprintf(fname,"%s/C%06d.fasta",dirname,k);
	T = fopen(fname,"w");
	k++;
      }
      for (j=i; j>=0; j = tree[j].next) 
	 fasta(T,j);
      if (count>=d) fclose(T);
    }
  }
  fclose(S0);
}




//-------------  auxiliary routines for common word matches

#ifndef NOINLINE
inline
#endif




void pseqi(int i) {  // depbugging code only
  // prints out the sequence (since sequences are stored in
  // compressed form, can't just do a printf!
  int k,curr;
  int  w=(SEQELTWIDTH/2);
  
  curr = i<0 ? -i : i;
  fprintf(outf,">%s\n",seqID[curr].id);
  if (i<0) {
    fprintf(outf,"NNNNNNNNNNNNNNNN\n");
    return;
  }
  i=curr;
  for(k=seqInfo[i].len-1; k>=0; k--) {
    fprintf(outf,"%c",cod2ch((~seq[i][k/w])>>((k%w)*2)&3));
    if ((seqInfo[i].len-k)%60==0) fprintf(outf,"\n");
  }
  if ((k%60)!=0) fprintf(outf,"\n");  
}


void pseq(int i) {  // depbugging code only
  // prints out the RC of the sequence
  int k,curr;
  int  w=(SEQELTWIDTH/2);

  curr = i<0 ? -i : i;
  fprintf(outf,">%s\n",seqID[curr].id);
  if (i<0) {
    fprintf(outf,"NNNNNNNNNNNNNNNN\n");
    return;
  }
  i=curr;
  for(k=0; k<seqInfo[i].len; k++) {
    fprintf(outf,"%c",cod2ch(seq[i][k/w]>>((k%w)*2)&3));
    if ((k+1)%60==0) fprintf(outf,"\n");
  }
  if ((k%60)!=0) fprintf(outf,"\n");  
}


void sample_heuristic(WorkPtr work, int s1, int s2, int *posmat, int *rcmat) {

  int i, limit ;
  uint16_t w0;
  unsigned short neg, pos;


  pos=neg=0;
  if (word_threshold < 5 ) return;
  NUM_num_matches++;

  limit = seqInfo[s2].len/BASESPERELT-1;
  for(i=1; i < limit ;i=i+2) {
      if (work->tableP[seq[s2][i]]) pos ++; 
      if (work->tableR[seq[s2][i]]) neg ++;
  }

  *posmat=pos;
  *rcmat =neg;
  //DBGprintf("Pos is %d\n",pos);
  if ((pos<=2) && (neg<= 2)) return;

  for(i=0; i < limit;i=i+2) {
    w0  =  seq[s2][i];
    if (work->tableP[w0]) pos++;
    if (work->tableR[w0]) neg++;
  }
  *posmat=pos;
  *rcmat =neg;
}


int tv_heuristic_pos(WorkPtr work, int s1, int s2) {
  // s1, s2 are the sequences
  int i,j;            // current character
  int np=0,npmax=0;   // how many common words
  uint16_t w1;

  NUM_full_heuristic++;
  w1=0;
  np=0;
  for(i=0; i<SAMPLE_WORD_LEN-1; i++) {
    w1 = GETWORD(w1, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
  }
  j=-window_len;
  for(i=SAMPLE_WORD_LEN-1; i < seqInfo[s2].len;i=i+1) {
    w1  =  GETWORD(w1, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
    work->tableQ[i]=work->tableP[w1];
    if (work->tableP[w1]) np++; 
    if ((j>0)  && work->tableQ[j]) np--;
    if (np>npmax) npmax=np;
    j++;
  }
  return npmax;
}


int tv_heuristic_rc(WorkPtr work, int s1, int s2) {
  // s1, s2 are the sequences
  // posmat is the number of common words between s1 and s2
  // rcmact is the number of common words between s1 and rc(s2)
  int i,j;            // current character
  int np=0,npmax=0;   // how many common words
  uint16_t w1;

  NUM_full_heuristic++;
  w1=0;
  np=0;
  for(i=0; i<SAMPLE_WORD_LEN-1; i++) {
    w1 = GETWORDRC(w1, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
  }
  j=-window_len;
  for(i=SAMPLE_WORD_LEN-1; i < seqInfo[s2].len;i=i+1) {
    w1  =  GETWORDRC(w1, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
    work->tableQ[i]=work->tableP[w1];
    if (work->tableP[w1]) np++; 
    if ((j>0)  && work->tableQ[j]) np--;
    if (np>npmax) npmax=np;
    j++;
  }
  return npmax;
}



// kept for old times' sake
void tv_heuristic(WorkPtr work, int s1, int s2, int *posmat, int *rcmat) {
  // s1, s2 are the sequences
  // posmat is the number of common words between s1 and s2
  // rcmact is the number of common words between s1 and rc(s2)
  int i;            // current character
  int k,prev3, prev2;
  int np=0,nrc=0;   // how many common words
  uint16_t w1,w1n;

  *posmat=*rcmat=0;
  prev3=prev2= (-2-GAP);

  NUM_full_heuristic++;
  k=w1=w1n=0;
  np=nrc=0;
  for(i=0; i<SAMPLE_WORD_LEN-1; i++) {
    w1 = GETWORD(w1, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
    w1n= GETWORDRC(w1n, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);

  }
  for(i=SAMPLE_WORD_LEN-1; i < seqInfo[s2].len;i=i+1) {
    w1  =  GETWORD(w1, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
    w1n =  GETWORDRC(w1n, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
    //DBGprintf("i %d: %d-%d %s %s\n",i, w1n,rc_big[w1],word2str(w1,8),word2str(w1n,8));
    if (work->tableP[w1] && i-prev2 >= GAP) {
      if (i-prev2 <= window_len)
	np++;
      else
	if (np < word_threshold) np=1;
      prev2=i;
      //continue;
    }
    if (work->tableP[rc_big[w1]] && i-prev3 >= GAP) {
      if (i-prev3 <= window_len) 
	nrc++;
      else
	if (nrc < word_threshold) nrc = 1;
      prev3 = i;
    }
  }
  *posmat = np;
  *rcmat = nrc;
}



void stat_num_matches(WorkPtr work, int s1, int s2, int *posmat, int *rcmat) {
  // s1, s2 are the sequences
  // posmat is the number of common words between s1 and s2
  // rcmact is the number of common words between s1 and rc(s2)
  int i,k;            // current character
  int prev3, prev2;
  int np=0,nrc=0, pmax=0, rmax=0;   // how many common words
  wordElt w2, w3;

  *posmat=*rcmat=0;

  prev3=prev2= (-2-GAP);
  w2=0;
  NUM_num_matches++;

  // we first do a sampling to see whether it's even worth 
  // checking the number of common words
  // NB: we are conservative, so if the user chooses a small
  // word_threshold, we do a comprehensive search
  if (word_threshold >= 0 ) {
    for(i=0; i < seqInfo[s2].len/BASESPERELT;i=i+1) {
      w2  =  seq[s2][i]&BIG_WORD_MASK;
      w3  =  (seq[s2][i]>>12)&BIG_WORD_MASK;
      if (work->tableP[w2]) np++;
      if (work->tableP[rc_big[w3]]) nrc++;
    }
  }
  // now if it looks plausible, we do an in-depth check
  sp=np;  // note horrible kludge sp, sr are global
  sr=nrc;
  if ( (word_threshold < 5) || (np >= 0) || (nrc >= 0) ) {
    k=w2=0;
    for(i=0; i<SAMPLE_WORD_LEN-1; i++) {
      w2 = GETWORD(w2, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
      k++;
    }
    np=nrc=0;
    for(i=SAMPLE_WORD_LEN-1; i < seqInfo[s2].len-(SAMPLE_WORD_LEN-word_len);i=i+1) {
      k++;
      w2  =  GETWORD(w2, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
      if (work->tableP[w2] && i-prev2 >= GAP) {
        if (i-prev2 <= window_len)
  	  np++;
	else
	  np=1;
	if (np>pmax) pmax=np;
	prev2=i;
	continue;
      }
      if (work->tableP[rc_big[w2]] && i-prev3>=GAP) {
	if (i-prev3 <= window_len) 
	  nrc++;
	else
	  nrc = 1;
        if(nrc>rmax) rmax=nrc;
	prev3 = i;
	continue;
      }
    }
  }
  *posmat = pmax;
  *rcmat = rmax;
}



//----------------------------------------------------------



void initialise(WorkPtr work, FILE * edfile) {
  int i,j,wlen;
  // Initialises all the main global variables
  // setting constants-----------------------

  wlen = word_len;
  tr[0]=0;   // to represent the bases. Bits 2+3 of the ASCII
  tr[1]=1;   // code give an encoding  A,C,T,G, but I want
  tr[2]=3;   // A,C,G,T so that ~ gives me the complement
  tr[3]=2;
  word_tsize = (1 << (2*wlen));  // word_tsize = 4^{word_len}
  word_mask  = word_tsize-1;  
  word_shift = 2*(wlen-1);
  for(j=0; j<num_seqs; j++) {
      tree[j].cluster  = j;
      tree[j].last     = j;
      tree[j].next     = -1;
      tree[j].match    = -1;
      tree[j].orient   = 1;
      tree[j].rank     = 0;
  }
  //----------memory allocation
  for(i=0; i<num_threads; i++) {
    work[i].workflag=0;
    work[i].thread_num = i;
    work[i].tableP = (short *) calloc(BIG_WORD_TSIZE+32, sizeof(short));
    work[i].tableR = (short *) calloc(BIG_WORD_TSIZE+32, sizeof(short));
    bzero(work[i].tableP,(BIG_WORD_TSIZE+32)*sizeof(short));
    bzero(work[i].tableR,(BIG_WORD_TSIZE+32)*sizeof(short));
    work[i].word1  = (wordElt *) calloc(MAX_SEQ_LEN, sizeof(wordElt));
    work[i].word2  = (wordElt *) calloc(MAX_SEQ_LEN, sizeof(wordElt));
  }

  if (dist == ed)   edinit(work,edfile);
  else  d2init(work);


  make_rc(rc_big, SAMPLE_WORD_LEN, BIG_WORD_TSIZE);
  // now allocate word_tsize integers: default an array of 4096 ints
  rc = (SeqElt *) calloc(word_tsize, sizeof(int));
  make_rc(rc,     wlen,     word_tsize);
}


int b2e(int numbytes) {
  // returns number of seqElt
  return (numbytes+4)/BASESPERELT+1;
}

int count_seqs(char * fname, int *datasize) {
  int i, count=0, total=0, curr=0,m;
  char dummy[MAX_SEQ_LEN];
  FILE * finp;

  finp = fopen(fname,"r");
  if (finp == 0)  {
    perror(fname);
    exit(5);
  }

  /* Store the sequence as an array of words */
  for(i=0; count<num_seqs || num_seqs==0; ) {
    m=fscanf(finp,"%[>]",dummy);
    if (m==EOF) break;
    if (m==1) { 
      count++;
      if (curr>0) {
        // adjust so we only allocate in block sizes of 16 byes
	//printf("%d bytes %d blocks\n",curr,b2e(curr));
	total = b2e(curr)+total;
	curr=0;
      }	  
      m=fscanf(finp,"%[] \t0-9@#$^&*()_=+;:'\"|\%,<.>/?`~A-Za-z{}[-]\n",
	       dummy);
    } else {
      m=fscanf(finp,"%[] \t0-9@#$^&*()_=+;:'\"|\%,<.>/?`~A-Za-z{}[-]\n",
	       dummy);
      curr = curr+strlen(dummy);
    }

  }
  total = b2e(curr)+total;
  *datasize = total;
  return count;
}




void output_rc(FILE * fout, char *line) {
  int i,l;
  char x;
  l = strlen(line);
  for(i=l-1; i>=0; i--) {
    x=line[i];
    switch(line[i]) {
    case 'A': x='T'; break;
    case 'a': x='t'; break;
    case 'C': x='G'; break;
    case 'c': x='g'; break;
    case 'G': x='C'; break;
    case 'g': x='c'; break;
    case 'T': x='A'; break;
    case 't': x='a'; break;
    }
    fprintf(fout,"%c",x);
    if ((l-i)%60==0) fprintf(fout,"\n");
  }
  if ((l-i-1)%60 != 0) fprintf(fout,"\n");
}

void output_fwd(FILE * fout, char *line) {
  int i,l;
  char x;
  l = strlen(line);
  for(i=0; i<l; i++) {
    fprintf(fout,"%c",line[i]);
    if ((i+1)%60==0) fprintf(fout,"\n");
  }
  if (i%60 != 0) fprintf(fout,"\n");
}


inline
void chomp(char * line) {
  int n;
  n=strlen(line)-1;
  if (line[n]='\n')
    line[n]=(char) 0;
}

void mirror_sequences(char *infname) {
  char empty[32];
  int i=0;
  for(i=0; i<num_seqs; i++) 
    pseq(IGNORE_SEQ(i)?-i:i);
  for(i=0; i<num_seqs; i++) 
    pseqi(IGNORE_SEQ(i)?-i:i);
}

void read_sequences(FILE *finp, int c1, int c2) {
  /* finp -- file we read from
     sequences are read into seq[c1],...,seq[c2-1] */
  int i, j, m,k,shift;
  char next_seq[MAX_SEQ_LEN],dummy[MAX_SEQ_LEN],clone[256];
  char tmp[MAX_SEQ_LEN],id_seq[4096];
  SeqPtr   curr;

  /* Store the sequence as an array of words */

  // where in data array do we start
  if   (c1==0) curr =  data;
  else
    curr = seq[c1-1]+(seqInfo[c1-1].len/BASESPERELT+4);
  for(i=c1; i<c2; i++) {
    seq[i] = curr;
    m=fscanf(finp,"%[>]",dummy);
    if (m != 1) {
      printf("Error in reading sequences: number %d %d of %d\n", i, m,num_seqs);
      fscanf(finp,"%s",dummy);
      printf("Label should start with a '>', found %s\n", dummy);
      printf("Last sequence successfully read in was: %s\n", id_seq);
      assert(0);
    }
    fscanf(finp,"%*[ \t]");
    m=fscanf(finp,"%s",id_seq);
    strcpy(clone,"");
    if (m==0) sprintf(id_seq,"%d",i);
    else {
      while (m=fscanf(finp,"%[\n]", tmp)==0) {
	m=fscanf(finp,"%[\t ]",dummy);
	m=fscanf(finp,"%[]\r0-9@#$^&*()_=+;:'\"|\%,<.>/?`~A-Za-z{}[-]",dummy);
	if (strcmp(dummy,"clone")==0) {
	  m = fscanf(finp," %[]@#$^&\\*()_=+;:'\"|\%,<.>/?`~0-9A-Za-z{}[-]",clone);
	  if (m != 1)  {
	    printf("Error in reading : %s # %d of %d\n",id_seq,i,num_seqs);
	    printf("clone keyword but no clone id\n");
	    assert(0);
	  }
	}
      }
    }
    m = fscanf(finp,"\n");
    strcpy(next_seq,"");
    while (fscanf(finp,"%[-A-Za-z0-9]\n",tmp)==1)  {
      strcat(next_seq,tmp);
    }
    fscanf(finp,"%*[.\n]");
    k = 1+strlen(next_seq)*2/SEQELTWIDTH;
    seqInfo[i].len = strlen(next_seq);
    tree[i].cluster  = i;
    tree[i].last     = i;
    tree[i].next     = -1;
    tree[i].match    = -1;
    tree[i].orient   = 1;
    tree[i].rank     = 0;
    if (IGNORE_SEQ(i)) continue;
#ifndef NOAUXINFO    
    seqID[i].clone = (string)malloc(strlen(clone)+1);
    strcpy(seqID[i].clone,clone);
    seqID[i].id       = (string) malloc(strlen(id_seq)+1);
    strcpy(seqID[i].id,id_seq);
    //printf("%s seq[%d] len=%d\n",seqID[i].id,i,seqInfo[i].len);
#endif
    shift=k=0;
    for(j=0; j< seqInfo[i].len; j++) {
      seq[i][k] = (seq[i][k]>>2) | (code(next_seq[j])<<BIG_SHIFT) ;
      if (shift == SEQELTWIDTH-2) {
	k++;
	shift= -2;
      }
      shift = shift+2;
    }
    for(j=seqInfo[i].len; j%BASESPERELT !=0 ; j++)
      curr[k] = curr[k]>>2;
    curr= curr+b2e(seqInfo[i].len);
   }
}


void init_sequences(FILE *finp, int c1, int c2) {

  int i;
  for(i=c1; i<c2; i++) {
    seqInfo[i].len =-1;
    seq[i] = (SeqPtr)  0;
    tree[i].cluster  = i;
    tree[i].last     = i;
    tree[i].next     = -1;
    tree[i].match    = -1;
    tree[i].orient   = 1;
    tree[i].rank     = 0;
    if (IGNORE_SEQ(i)) continue;
  }
}



void show_EXT(FILE * outf) {
  int i;
  fprintf(outf,"   Index      Cluster     Link   Orientation   Witness\n");
  for (reindex_value=0; i<num_seqs; i++) 
    fprintf(outf,"%10d %10d %10d %10d %10d\n",
	    i,find_parent(i),tree[i].next,tree[i].orient,tree[i].match);
}




void create_word_lists(WorkPtr work, int s1, int s2, int rcflag ) {
  // refers to global variables word1 & word2
  // procedure uncompresses the two sequences
  int w1,w2;
  int i;
  //init
  w1=w2=0;
  for(i=0; i<word_len-1; i++) {
    w1 = GETWORD(w1, seq[s1], i,word_shift, word_mask);
    w2 = GETWORD(w2, seq[s2], i, word_shift, word_mask);
  }
  //now fill in the words 1 by 1
  for(i=word_len-1; i < seqInfo[s1].len; i++) {
    w1 = GETWORD(w1, seq[s1], i,word_shift, word_mask);
    //DBGprintf(">w1 %d %d %s %s\n",i, w1, word2str(w1,6),word2str(rc[w1],6));
    work->word1[i-word_len+1]=w1;
  }
  for(i=word_len-1; i < seqInfo[s2].len; i++) {
    w2 =  GETWORD(w2, seq[s2], i, word_shift, word_mask);
    work->word2 [i-word_len+1]= rcflag ? rc[w2] : w2;
  }
}



#define BIG_WORD_LEN 8


void get_bounds(WorkPtr work, int s1, int s2, int * left, int * right) {
  // Once we do decide to compare s1 and s2 it would still be
  // wasteful to compare the whole of s1 with the whole of s2.
  // Thus we find the regions of s1 which might match with something
  // with s2: left and right are set to these values
  SeqElt w, w2;
  int i;
  //init
  *left = -1;
  *right= -1;
  w=w2=0;
  if (word_threshold == 0) {
    *left=0;
    *right=seqInfo[s1].len;
    return;
  }
  for(i=0; i<BIG_WORD_LEN-1; i++) {
    w = GETWORD(w, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
  }
  for(i=BIG_WORD_LEN-1; i < seqInfo[s2].len; i++) {
    w = GETWORD(w, seq[s2], i,BIG_SHIFT, BIG_WORD_MASK);
    work->word_seenQ[i-BIG_WORD_LEN+1]=w;
    work->tableQ[w]=1;
  }
  w=w2=0;
  for(i=0; i<BIG_WORD_LEN-1; i++) {
    w = GETWORD(w, seq[s1], i,BIG_WORD_LEN, BIG_WORD_MASK);
    w2= GETWORDRC(w2, seq[s1], i,BIG_SHIFT, BIG_WORD_MASK);
  }
  for(i=BIG_WORD_LEN-1; i < seqInfo[s1].len; i++) {
    w =  GETWORD(w, seq[s1], i, BIG_SHIFT, BIG_WORD_MASK);
    w2= GETWORDRC(w2, seq[s1], i,BIG_SHIFT, BIG_WORD_MASK);
    if (work->tableQ[w] || work->tableQ[w2]) {
      *right = i;
      if (*left < 0) *left = i;
    }
  }
  for(i=0; i< seqInfo[s2].len; i++) 
    work->tableQ[work->word_seenQ[i]]=0;
}


int get_next_bounds(WorkPtr work,int s1, int s2, int * left, int * right) {
  // Once we do decide to compare s1 and s2 it would still be
  // wasteful to compare the whole of s1 with the whole of s2.
  // Thus we find the regions of s1 which might match with something
  // with s2: left and right are set to these values
  SeqElt w, w2;
  int i;
  //init
  *left = -1;
  *right= -1;
  w=w2=0;
  if (word_threshold == 0) {
    *left=0;
    *right=seqInfo[s1].len;
    return 1;
  }
  for(i=0; i<BIG_WORD_LEN-1; i++) {
    w = GETWORD(w, seq[s2], i,BIG_WORD_LEN, BIG_WORD_MASK);
  }
  for(i=BIG_WORD_LEN-1; i < seqInfo[s2].len; i++) {
    w = GETWORD(w, seq[s2], i,word_shift, BIG_WORD_MASK);
    work->word_seenQ[i-BIG_WORD_LEN+1]=w;
    work->tableQ[w]=1;
  }
  w=w2=0;
  for(i=0; i<BIG_WORD_LEN-1; i++) {
    w = GETWORD(w, seq[s1], i,BIG_WORD_LEN, BIG_WORD_MASK);
    w2= ((((~w)&3)<< BIG_SHIFT)) | (w2 >> 2);
  }
  for(i=BIG_WORD_LEN-1; i < seqInfo[s1].len; i++) {
    w =  GETWORD(w, seq[s1], i, word_shift, BIG_WORD_MASK);
    w2= ((((~w)&3)<< BIG_SHIFT)) | (w2 >> 2);
    if (work->tableQ[w] || work->tableQ[w2]) {
      *right = i;
      if (*left < 0) *left = i;
    }
  }
  for(i=0; i< seqInfo[s2].len; i++) 
    work->tableQ[work->word_seenQ[i]]=0;
  return 1;
}

