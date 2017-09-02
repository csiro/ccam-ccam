#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
 
void getstacklimit_(int *limit)
{
    struct rlimit rlim;
    int ierr;
    ierr =  getrlimit(RLIMIT_STACK, &rlim);
    /* printf("Stack limit %d %d\n", rlim.rlim_max, rlim.rlim_cur); */
    *limit = rlim.rlim_cur;
}

void setstacklimit_(int *limit)
{
    struct rlimit rlim;
    int ierr;
    rlim.rlim_max = *limit;
    rlim.rlim_cur = *limit;
    ierr =  setrlimit(RLIMIT_STACK, &rlim);
    /* printf("Set stack limit %d \n", ierr); */
}
