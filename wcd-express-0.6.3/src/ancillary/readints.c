
#include <getopt.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>


static char * allowed_options = "f:t:i:";


unsigned long long int from, until;


void process_options(int argc, char * argv[]) {
  int opt, m;
  char outfname[255], resp;

  while ((opt = getopt(argc,argv,allowed_options)) != -1) {
    switch(opt) {
    case 'f':
      sscanf(optarg,"%lld",&from);
      break;
    case 't':
      sscanf(optarg,"%lld",&until);
      break;
    case 'i':
      sscanf(optarg,"%lld",&until);
      from = until;
      break;
    }
  }


}


int main(int argc, char *argv[]) {

  unsigned long long int i, num;
  int fd,n;

  from=until=0;
   
  process_options(argc,argv);
  fd = open(argv[optind],0,0);

  if (fd<0) {
    perror(argv[optind]);
    exit(1);
  }

  i =0;
  do {
    n=read(fd,&num,8);
    if ((n<8) || (i>until)&&(until>0))  break;
    if (i>=from)
      printf("%10lld %10lld\n",i,num);
    i++;
  } while (1);
  close(fd);

}
