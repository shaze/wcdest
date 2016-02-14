bin_PROGRAMS = wcd
wcd_SOURCES = wcd.c d2.c ed.c common.c wcd.h d2.h ed.h common.h suffixcluster.c suffixcluster.h mpistuff.h mpistuff.c pthreadstuff.c pthreadstuff.h auxcluster.h auxcluster.c
wcd_CFLAGS = 
wcd_LDFLAGS = 
if MPI
wcd_CFLAGS += -DMPI
endif
if PTHREADS
wcd_CFLAGS += -DPTHREADS
wcd_LDFLAGS += -lpthread
endif
if !COMPILER
wcd_CFLAGS += -DNOINLINE 
endif
if !LONGOPT
wcd_CFLAGS += -DNOLONGOPT
endif
if !LONGOPT2
wcd_CFLAGS += -DNOLONGOPT
endif
