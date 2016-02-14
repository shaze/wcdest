


#include "common.h"
#include "strings.h"
#include "assert.h"
#include "string.h"

extern SeqInfoPtr seqInfo;

extern int num_seqs;
extern int skip;
extern FILE *outf;
extern int ignore_active;


void process_split(char * fname, char * prefix) {
  FILE * finp, * fout;
  int    res, cnum, cluster, j;
  char  outname[255];
  int   *curr_cluster;

  curr_cluster = calloc(num_seqs,sizeof(int));

  cluster = 0;
  finp = fopen(fname,"r");
  if (finp == NULL)  {
    perror(fname);
    exit(1);
  }
  while (fscanf(finp,"%d",&cnum)==1) {
    j=0;
    do {
      curr_cluster[j]=cnum;
      j++;
    } while(fscanf(finp, "%d", &cnum)==1);
    if (j>=skip) { // NB: we are abusing the name skip
                   // in this context it's just the threshold for 
                   // cluster size
      sprintf(outname,"%s%05d.fa",prefix,cluster);
      cluster++;
      fout = fopen(outname,"w");
      for(cnum=0; cnum<j; cnum++)
	fasta(fout, curr_cluster[cnum]);
      fclose(fout);
    }
    res = fscanf(finp,"%[.]",outname);
    if (res != 1) {
      strcpy(outname,"PREMATURE END OF LINE");
      fscanf(finp,"%s",outname);
      fprintf(outf,"Error in split file <%s>: expected '.' found <%s>\n",
	      fname,outname);
      assert(0);
    }
  }
}


void process_constraints(char * fname, int base) {
  FILE * finp;
  int    res, cnum, flag_val,i,invert=0;
  char   command[255];
  int plusrcflag=0;
  finp = fopen(fname,"r");



  if (finp == NULL)  {
    perror(fname);
    exit(1);
  }
  while (fscanf(finp,"%s",command)==1) {

    if (strcasecmp(command,"fix")==0) flag_val = FIX; else
      if (strcasecmp(command,"cluster-only")==0) {
        flag_val = IGNORE;
	ignore_active = 1;
      } else 
	if (strcasecmp(command,"cluster-others")==0) {
	  flag_val = IGNORE;
	  ignore_active = 0;
	} else
	  if (strcasecmp(command,"reset-others")==0) {
	    for(i=0; i<num_seqs; i++) SET_FLAG(RESET,i);
	    ignore_active=1;
	    invert=1;
	  } else
	    if (strcasecmp(command,"reset")==0) {
	      flag_val = RESET;
	      ignore_active = 1;
	    } else 
	      if (strcasecmp(command,"cluster-others.")==0) {
		return; //empty -- for automated pipelines which
                        // produce filters that may be empty
	      } else {
		fprintf(outf,"Illegal value in file %s. Unrecognised constraint type <%s>\n",
			fname, command);
		assert(0);
	      }
    while(fscanf(finp, "%d", &cnum)==1) {
      if (invert)
	seqInfo[cnum+base].flag=0;
      else {
	SET_FLAG(flag_val, cnum+base);
        if (plusrcflag) // We have file followed by RC
	  SET_FLAG(flag_val, cnum+base+num_seqs/2);
      }
    }
    
    res = fscanf(finp,"%[.]",command);
    if (res != 1) {
      strcpy(command,"PREMATURE END OF LINE");
      fscanf(finp,"%s",command);
      fprintf(outf,"Error in constraint file <%s>: expected '.' found <%s>\n",
	      fname,command);
      assert(0);
    }
    
  }

}

