wcd is an open source tool for clustering expressed sequence tags. It implements matching using either d2 or edit distance, which provides sensitive matching. Heuristics speed up the clustering process. wcd parallelises well under MPI and also provides Pthreads clustering.
wcd has a very compact memory representation which means that large files can be clustered on a single processor (and it can run quite happily without MPI or Pthreads headers or libraries being available). Of course, if you have a 100 CPU cluster -- it'll go much faster!