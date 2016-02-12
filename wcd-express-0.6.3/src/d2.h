

#define __D2HINCL

#ifndef __COMMONINCL
#include "common.h"
#endif

int  d2(WorkPtr work, int s1, int s2, int rcflag);

int  d2pair(WorkPtr work, int s1, int s2, int rcflag);

int  endPoints(WorkPtr work, int s1, int s2, int rcflag);

void d2init(WorkPtr threadList);


