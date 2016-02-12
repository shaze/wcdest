
/* 	$Id: ed.c,v 1.1 2007/10/19 14:40:33 scott Exp scott $	 */

#ifndef lint
static char vcid[] = "$Id: ed.c,v 1.1 2007/10/19 14:40:33 scott Exp scott $";
     #endif /* lint */

/*

Copyright (C) Scott Hazelhurst     2003, 2004
              School of Computer Science, 
              University of the Witwatersrand,
	      Johannesburg
	      Private Bag 3, 2050 Wits
              South Africa
              scott@cs.wits.ac.za

     This program is free software; you can redistribute it and/or
     modify it under the terms of the GNU General Public Licence
     as published by the Free Software Foundation; either version 2
     of the Licence, or (at your option) any later version.
     
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public Licence for more details.

      ftp://ftp.cs.wits.ac.za/pub/research/software/wcd_distrib.tar.gz

*/

#include <string.h>
#include "ed.h"
#include "wcd.h"

typedef struct dmatEntry  {
  short e,f, g,m;
} dmatEntry;

typedef dmatEntry * dmat1Arr;



static int cost[5][5];


 // imported from common.c

extern int NUM_dfn;
extern int word_len;
extern int window_len;
extern int num_threads;

// impored from wcd.c

extern SeqInfoPtr seqInfo;

extern int theta; // default value

// default penalties for distance function
// edinit might change these
static int EXTEND=1; // cost when gap extends
static int OPEN  =1; // cost when gap opens
static int MISMATCH=2;
static int MATCH=-1;

#define MIN3(a,b,c) a < b ? (a<c?a:c) : (b<c?b:c)

int  ed(WorkPtr work, int s1, int s2, int rcflag) {
  // if rcflag==1 then we study the rc(s2) rather than s2
  int row, col, l1, l2, r1, r2, curr,prev;

  int num_words_1, num_words_2;
  int score=0;
  void *x, *y;
  dmat1Arr dmat[2];
  int the_min,min_min=0;


  dmat[0] = work->fn_data;
  dmat[1] = work->fn_data+MAX_SEQ_LEN;

  //memset(dmat[0],0,2*sizeof(dmatEntry)*MAX_SEQ_LEN);

  // rather than comparing the entire ESTs, we narrow down
  // the regions on the ESTs where the match could be
  get_bounds(work,s1, s2, &l1, &r1);  // bounds for s1
  get_bounds(work,s2, s1, &l2, &r2);  // bounds for s2
  if (r1-l1< 4)return (10*window_len);// don't bother if not there
  l1 = MAX(0,l1-window_len);     // region could star 
  l2 = MAX(0,l2-window_len);
  r1 = MIN(r1+window_len,seqInfo[s1].len);
  r2 = MIN(r2+window_len,seqInfo[s2].len);
  NUM_dfn++;
  num_words_1 = r1-word_len+1;
  num_words_2 = r2-word_len+1;

  create_word_lists(work,s1, s2, rcflag);
  prev=0;
  curr=1;
  for(row=0; row<=r1; row++) dmat[0][row].f = OPEN;
  for(col=l2; col<r2 && min_min  > theta ; col++) {
    dmat[curr][l1].f=dmat[curr][l1].m=0;
    dmat[curr][l1].e=OPEN;
    for (row=l1+1; row<=r1 && min_min > theta; row++) {
       the_min=0;
       dmat[curr][row].f = 
          MIN3(dmat[prev][row].e+EXTEND+OPEN,
               dmat[prev][row].f+EXTEND,
               dmat[prev][row].g+OPEN+EXTEND);
       if (dmat[curr][row].f < the_min) the_min = dmat[curr][row].f;
       dmat[curr][row].g = dmat[prev][row-1].m+
	 cost[work->word1[row-1]][work->word2[col]];
       if (dmat[curr][row].g < the_min) the_min = dmat[curr][row].g;
       dmat[curr][row].e = 
          MIN3(dmat[curr][row-1].e+EXTEND,
               dmat[curr][row-1].f+EXTEND+OPEN,
               dmat[curr][row-1].g+OPEN+EXTEND);
       if (dmat[curr][row].e < the_min) the_min = dmat[curr][row].e;
       dmat[curr][row].m=the_min;
       if (the_min<min_min) min_min=the_min;
    }
    prev = curr;
    curr = 1-curr;
  }

  return min_min;
}


int  edpair(WorkPtr work, int s1, int s2, int rcflag) {
  // if rcflag==1 then we study the rc(s2) rather than s2
  int row, col, l1, l2, r1, r2, curr,prev;

  int num_words_1, num_words_2;
  int score=0;
  void *x, *y;
  dmat1Arr dmat[2];
  int the_min,min_min=0;


  dmat[0] = work->fn_data;
  dmat[1] = work->fn_data+MAX_SEQ_LEN;

  //memset(dmat[0],0,2*sizeof(dmatEntry)*MAX_SEQ_LEN);

  l1 = 0;     // region could star 
  l2 = 0;
  r1 = seqInfo[s1].len;
  r2 = seqInfo[s2].len;
  NUM_dfn++;
  num_words_1 = r1-word_len+1;
  num_words_2 = r2-word_len+1;

  create_word_lists(work,s1, s2, rcflag);
  prev=0;
  curr=1;
  for(row=0; row<=r1; row++) dmat[0][row].f = OPEN;
  for(col=l2; col<r2 && min_min  > theta ; col++) {
    dmat[curr][l1].f=dmat[curr][l1].m=0;
    dmat[curr][l1].e=OPEN;
    for (row=l1+1; row<=r1 && min_min > theta; row++) {
       the_min=0;
       dmat[curr][row].f = 
          MIN3(dmat[prev][row].e+EXTEND+OPEN,
               dmat[prev][row].f+EXTEND,
               dmat[prev][row].g+OPEN+EXTEND);
       if (dmat[curr][row].f < the_min) the_min = dmat[curr][row].f;
       dmat[curr][row].g = dmat[prev][row-1].m+
	 cost[work->word1[row-1]][work->word2[col]];
       if (dmat[curr][row].g < the_min) the_min = dmat[curr][row].g;
       dmat[curr][row].e = 
          MIN3(dmat[curr][row-1].e+EXTEND,
               dmat[curr][row-1].f+EXTEND+OPEN,
               dmat[curr][row-1].g+OPEN+EXTEND);
       if (dmat[curr][row].e < the_min) the_min = dmat[curr][row].e;
       dmat[curr][row].m=the_min;
       if (the_min<min_min) min_min=the_min;
    }
    prev = curr;
    curr = 1-curr;
  }

  return min_min;
}



void edinit(WorkPtr threadList, FILE *edfile) {
  int i, j;
  if (edfile == NULL) 
    for(i=0; i<4; i++)
      for(j=0;j<4;j++)
	cost[i][j] = i==j ? MATCH : MISMATCH;
  else {
    for(i=0; i<4; i++) 
      fscanf(edfile,"%d %d %d %d\n",
	     &cost[i][0],&cost[i][1],&cost[i][2],&cost[i][3]);
    fscanf(edfile,"%d %d",&OPEN,&EXTEND); 
  }
  for(i=0; i<num_threads; i++)  {
    threadList[i].fn_data = calloc(MAX_SEQ_LEN*2,sizeof(dmatEntry)) ;
  }
}

