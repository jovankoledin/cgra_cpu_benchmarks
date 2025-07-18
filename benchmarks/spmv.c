/*
Based on algorithm described here:
http://www.cs.berkeley.edu/~mhoemmen/matrix-seminar/slides/UCB_sparse_tutorial_1.pdf
Copied from: https://github.com/breagen/MachSuite/blob/master/spmv/crs/spmv.c
*/

#include <stdlib.h>
#include <stdio.h>

// These constants valid for the IEEE 494 bus interconnect matrix
#define NNZ 1666
#define N 494

#define TYPE double

void spmv(TYPE val[NNZ], int32_t cols[NNZ], int32_t rowDelimiters[N + 1],
          TYPE vec[N], TYPE out[N]);

struct bench_args_t {
  TYPE val[NNZ];
  int32_t cols[NNZ];
  int32_t rowDelimiters[N+1];
  TYPE vec[N];
  TYPE out[N];
};
void spmv(TYPE val[NNZ], int32_t cols[NNZ], int32_t rowDelimiters[N+1], TYPE vec[N], TYPE out[N]){
    int i, j;
    TYPE sum, Si;

    spmv_1 : for(i = 0; i < N; i++){
        sum = 0; Si = 0;
        int tmp_begin = rowDelimiters[i];
        int tmp_end = rowDelimiters[i+1];
        spmv_2 : for (j = tmp_begin; j < tmp_end; j++){
            Si = val[j] * vec[cols[j]];
            sum = sum + Si;
        }
        out[i] = sum;
    }
}
