

// Routines that deal with pthreads stuff
// There is only one public routine


#ifdef PTHREADS
pthread_mutex_t get_task_mutex;
pthread_mutex_t write_union_mutex;
pthread_mutex_t invert_mutex;
pthread_mutex_t dump_mutex;
pthread_mutex_t find_parent_mutex;
#endif


void parallel_do_cluster_servant(WorkPtr work);
