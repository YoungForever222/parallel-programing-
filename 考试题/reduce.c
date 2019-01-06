#include <math.h>
#include <omp.h>
#define NUM_THREADS 4
#define N 10000
#define factor 1.0
int main(){
    double sum = 0.0,Pi = 0.0;
    int k;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel num_threads(NUM_THREADS)
    {
    printf("Hello! threadID=%d\n",omp_get_thread_num());
    }
    #pragma omp parallel for private(k) reduction(+:sum)
    for(k=0;k<N;k++){
        sum = sum + pow(-1,k) / (2.0*k+1.0);
    }
    Pi = 4.0 * sum;
    printf("Pi = %lf",Pi);
}