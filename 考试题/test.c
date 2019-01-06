#include <omp.h>
#define N 20
int main(){
    int i;
    float a[N],b[N];
    for(i=0;i<N;i++){
        a[i] = b[i] = i*1.0;
    }
    #pragma omp parallel shared(a,b) private(i)
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            for(i=1;i<16;i=i+3){
                a[i+3] = a[i]+b[i];
            }
            #pragma omp section
            for(i=2;i<16;i=i+3){
                a[i+3] = a[i]+b[i];
            }
            #pragma omp section
            for(i=3;i<16;i=i+3){
                a[i+3] = a[i]+b[i];
            }
        }
    }
    for(i=1;i<20;i++)
        printf("%f\t",a[i]);
}