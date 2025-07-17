/*
Copied from https://github.com/breagen/MachSuite/blob/master/fft/strided/
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FFT_SIZE 16
#define twoPI 6.28318530717959

void fft(double real[FFT_SIZE], double img[FFT_SIZE], double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2]);

struct bench_args_t {
        double real[FFT_SIZE];
        double img[FFT_SIZE];
        double real_twid[FFT_SIZE/2];
        double img_twid[FFT_SIZE/2];
};

void fft(double real[FFT_SIZE], double img[FFT_SIZE], double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2]){
    int even, odd, span, log, rootindex;
    double temp;
    log = 0;

    for(span=FFT_SIZE>>1; span; span>>=1, log++){
        for(odd=span; odd<FFT_SIZE; odd++){
            odd |= span;
            even = odd ^ span;

            temp = real[even] + real[odd];
            real[odd] = real[even] - real[odd];
            real[even] = temp;

            temp = img[even] + img[odd];
            img[odd] = img[even] - img[odd];
            img[even] = temp;

            rootindex = (even<<log) & (FFT_SIZE - 1);
            if(rootindex){
                temp = real_twid[rootindex] * real[odd] -
                    img_twid[rootindex]  * img[odd];
                img[odd] = real_twid[rootindex]*img[odd] +
                    img_twid[rootindex]*real[odd];
                real[odd] = temp;
            }
        }
    }
}

int main(){
    double data_x[FFT_SIZE];
    double data_y[FFT_SIZE];
    double img[FFT_SIZE];
    double real[FFT_SIZE];
    int i;

    //set up twiddles...
    double typed;
    int n, N;
    N = FFT_SIZE;

    //Pre-calc twiddles
    for(n=0; n<(N>>1); n++){
        typed = (double)(twoPI*n/N);
        real[n] = cos(typed);
        img[n] = (-1.0)*sin(typed);
    }

    //Init data
    for(i=0; i < FFT_SIZE; i++) {
        data_x[i] = (double)(i);
        data_y[i] = (double)(i);
    }

    fft(data_x, data_y, real, img);

    return 0;
} 
