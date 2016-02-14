
/* 	$Id: d2.c,v 0.2.5.1 2006/01/05 08:11:33 scott Exp scott $	 */

#ifndef lint
static char vcid[] = "$Id: d2.c,v 0.2.5.1 2006/01/05 08:11:33 scott Exp scott $";
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


#include "d2.h"
#include "wcd.h"
#include <stdlib.h>
#include <string.h>






// imported from common.c

//extern wordElt * word1, *word2;

// imported from wcd.c
extern int num_words_win;      // number of words in window
extern int theta; // default value
extern int window_len;
extern int word_len;
extern int word_tsize;  // 4^word_len
extern int word_mask;   // word_tszie1
extern int word_shift;  // 2*(word_len-1)
extern int NM_threshold;
extern int skip;
extern int NUM_num_matches;
extern int NUM_dfn;
extern int NUM_dfn_succ;
extern int NUM_full_heuristic;
extern int num_threads;

extern SeqInfoPtr seqInfo;
extern SeqPtr *seq;
extern int offset;

#ifndef NOINLINE
inline
#endif
void update_window_1(int * delta, wordElt * word1, 
                     int i, int *min_score, int *score) {
  // the word at int i falls out, and the word at 
  // i+window_len moves in
  wordElt w;
  ASSERT(i<MAX_SEQ_LEN);
  w = word1[i+num_words_win]&word_mask;
  *score = *score + 2*delta[w]+1;
  delta[w]++;
  w = word1[i]&word_mask;
  *score = *score - 2*delta[w]+1;
  delta[w]--;
  if (*score < *min_score) *min_score = *score;
}



#ifndef NOINLINE
inline
#endif
void update_left_window_2(int * delta,
			  wordElt * word2,
                          int j, int *min_score, int *score, int shift) {
  // the word at  j-num_words_win moves in
  // the word at  j moves out
  wordElt w;
  w = (word2[j-num_words_win])&word_mask;
  *score = *score - 2*delta[w]+1;
  delta[w]--;
  w = (word2[j])&word_mask;
  *score = *score + 2*delta[w]+1;
  delta[w]++;
  if (*score < *min_score) *min_score = *score;
}


#ifndef NOINLINE
inline
#endif
void update_right_window_2(int * delta,
			   wordElt * word2,
                           int j, int *min_score, int *score, int shift) {
  // the word at uint j falls out
  // j+num_words_win moves in
  wordElt w;
  w  = (word2[j+num_words_win]>>shift)&word_mask;
  *score = *score - 2*delta[w]+1;
  delta[w]--;
  w = (word2[j]>>shift) & word_mask;
  *score = *score + 2*delta[w]+1;
  delta[w]++;
  if (*score < *min_score) *min_score = *score;
}




int  endPoints(WorkPtr work, int s1, int s2, int rcflag) {
  // if rcflag==1 then we study the rc(s2) rather than s2
  int first_end_word, last_start_word; /* moved from global */
  int i, j, l1, l2, min_score, len;
  wordElt w;
  wordElt * word1, * word2;
  int *delta;
  int num_words_1, num_words_2;
  int score=0;


  len = MIN(seqInfo[s1].len,seqInfo[s2].len);
  word1 = work->word1;
  word2 = work->word2;
  delta = (int *) work->fn_data;

  num_words_win = len-word_len+1;
  memset(delta,0,sizeof(int)*word_tsize);
  memset(work->tableQ,0,seqInfo[s2].len);
  create_word_lists(work,s1, s2, rcflag);
  
  if (rcflag)
    l1=0;
  else
    l1=seqInfo[s1].len-len-1;
  l2=0;
  //printf("%d: l1=%d, l2=%d nwl=%d\n",rcflag,l1,l2,num_words_win);
  for(i=0; i<num_words_win; i++) {
    w = word1[l1+i]&word_mask;
    ASSERT(score>=0);
    score = score + 2*delta[w]+1;
    delta[w]++;
    w = word2[l2+i]&word_mask;
    score = score - 2*delta[w]+1;
    delta[w]--;
  }
  for(i=0; i<num_words_win; i++) {
    w = word1[l1+i]&word_mask;
    delta[w]=0;
    w = word2[l2+i]&word_mask;
    delta[w]=0;
  }
  min_score=score;
  score=0;
  if (rcflag)
    l1=seqInfo[s1].len-len-1;
  else
    l1=0;
  l2=seqInfo[s2].len-len-1;
  for(i=0; i<num_words_win; i++) {
    w = word1[l1+i]&word_mask;
    score = score + 2*delta[w]+1;
    delta[word1[l1+i]&word_mask]++;
    w = word2[l2+i]&word_mask;
    score = score - 2*delta[w]+1;
    delta[w]--;
    ASSERT(score>=0);
  }
  min_score=MIN(min_score,score);;
  return min_score;
}


int  d2(WorkPtr work, int s1, int s2, int rcflag) {
  // if rcflag==1 then we study the rc(s2) rather than s2
  int first_end_word, last_start_word; /* moved from global */
  int i, j, k, l1, l2, r1, r2, min_score;
  wordElt w;
  wordElt * word1, * word2;
  int *delta;
  int num_words_1, num_words_2;
  int score=0;
  int  shift;


  word1 = work->word1;
  word2 = work->word2;
  delta = (int *) work->fn_data;

  memset(delta,0,sizeof(int)*word_tsize);
  memset(work->tableQ,0,seqInfo[s2].len);

  num_words_win = window_len-word_len+1;

  // rather than comparing the entire ESTs, we narrow down
  // the regions on the ESTs where the match could be
  get_bounds(work,s1, s2, &l1, &r1);  // bounds for s1
  get_bounds(work,s2, s1, &l2, &r2);  // bounds for s2

  if (r1-l1< 4)return (10*theta);// don't bother if not there

  l1 = MAX(0,l1-window_len);     // region could star 
  l2 = MAX(0,l2-window_len);
  r1 = MIN(r1+window_len,seqInfo[s1].len);
  r2 = MIN(r2+window_len,seqInfo[s2].len);
  NUM_dfn++;
  num_words_1 = r1-word_len+1;
  num_words_2 = r2-word_len+1;

  create_word_lists(work,s1, s2, rcflag);

  // work out score for first windows-----------------
  for(i=0; i<num_words_win; i++) {
    score = score + 2*delta[word1[l1+i]&word_mask]+1;
    delta[word1[l1+i]&word_mask]++;
    w = word2[l2+i]&word_mask;
    score = score - 2*delta[w]+1;
    delta[w]--;
  }
  i=l1;
  // last_start_word: last word that can start a window in seq2
  last_start_word = num_words_2-num_words_win;
  // first_end_word:  first word that can end a window 
  first_end_word = l2+num_words_win-1;
  // now go through all the other pairs of windows
  min_score = score;
  //printf("s1=%d s2=%d s1.len=%d s2.len=%d l1=%d num_words_1=%d num_words_win=%d\n",
  //       s1,s2,seqInfo[s1].len,seqInfo[s2].len,l1,num_words_1,num_words_win);
  //return 100;
  do {
    // check every window in seq2 against curr window in seq1
    for(j=l2; j<last_start_word; j++) 
      update_right_window_2(delta,word2,j, &min_score, &score,0);
    if (min_score <= theta || i==num_words_1-num_words_win) break;
    //**if (i==num_words_1-num_words_win) break;
    k=0;
    do { // skip over a number of windows in seq1
      update_window_1(delta,word1,i, &min_score, & score);
      i++;
      if (score <= theta || i==num_words_1-num_words_win) break;
      //*if (i==num_words_1-num_words_win) break;
      k ++;
    } while(k < skip);
    // check every window in seq2 against curr window in seq1
    for(j=num_words_2-1; j > first_end_word ; j--)
      update_left_window_2(delta,word2,j, &min_score, &score,0);
    if (min_score <= theta ||i>=num_words_1-num_words_win) break;
    //if (i>=num_words_1-num_words_win) break;
    k=0;
    do { // skip over a number of windows in seq 1
      update_window_1(delta,word1,i, &min_score, & score);
      i++;
      if (score <= theta || i==num_words_1-num_words_win) break;
      //if (i==num_words_1-num_words_win) break;
      k ++;
    } while(k < skip);
  } while(1);
  if ((min_score>theta) && (min_score<0)) {
    memset(delta,0,sizeof(int)*word_tsize);
    score = endPoints(work,s1,s2,rcflag);
    // ** if (score<14) min_score=MIN(score,min_score);
  }
  return min_score;
}




//-------------------------------------------------------------
// Code to compare two sequences

//-------- code that prints d2 score of all pairs of windows between
//         two sequences. It's the same d2 scoer as above but
//         less efficient



int  d2pair(WorkPtr work, int s1, int s2, int rcflag) {
  // if rcflag==1 then we study the rc(s2) rather than s2
  int first_end_word, last_start_word; /* moved from global */
  int i, j, k;
  int num_words_1, num_words_2;
  int score=0,min_score;
  int * delta;
  wordElt * word1, *word2;

  word1 = work->word1;
  word2 = work->word2;
  delta = (int *) work->fn_data ;

  num_words_win = window_len-word_len+1;

  NUM_dfn++;
  memset(delta,0,sizeof(int)*word_tsize);

  num_words_1 = seqInfo[s1].len-word_len+1;
  num_words_2 = seqInfo[s2].len-word_len+1;

  create_word_lists(work,s1, s2, rcflag);

  // work out score for first windows-----------------
  for(i=0; i<num_words_win; i++) {
    score = score + 2*delta[word1[i]&word_mask]+1;
    delta[word1[i]&word_mask]++;
    score = score - 2*delta[word2[i]&word_mask]+1;
    delta[word2[i]&word_mask]--;
  }
  i=0;
  // last_start_word: last word that can start a window in seq2
  last_start_word = num_words_2-num_words_win;
  // first_end_word:  first word that can end a window 
  first_end_word = num_words_win-1;
  // now go through all the other pairs of windows
  min_score = score;
  do {
    // check every window in seq2 against curr window in seq1
    for(j=0; j<last_start_word; j++) {
      update_right_window_2(delta,word2,j, &min_score, &score,0);
      //if (score==0) printf("0 at %d %d\n",i,j);
    }
    if (i==num_words_1-num_words_win) return min_score;
    k=0;
    do { // skip over a number of windows in seq1
      update_window_1(delta,word1,i, &min_score, & score);
      //if (score==0) printf("0 at %d %d\n",i,j);
      i++;
      if (i==num_words_1-num_words_win) return min_score;
      k ++;
    } while(k < skip);
    // check every window in seq2 against curr window in seq1
    for(j=num_words_2-1; j > first_end_word; j--) {
      update_left_window_2(delta,word2,j, &min_score,&score,0);
      //if (score==0) printf("0 at %d %d\n",i,j);
    }
    if (i==num_words_1-num_words_win) return min_score;
    k=0;
    do { // skip over a number of windows in seq 1
      update_window_1(delta,word1,i, &min_score, & score);
      //if (score==0) printf("0 at %d %d\n",i,j);
      i++;
      if (i==num_words_1-num_words_win) return min_score;
      k ++;
    } while(k < skip);
  } while(1);
  assert(0);
  return 0;// should never get here
}


void d2init(WorkPtr threadList) {
  int i;
  for(i=0; i<num_threads; i++) {
     threadList[i].fn_data = calloc(word_tsize,sizeof(int));
  }
}


