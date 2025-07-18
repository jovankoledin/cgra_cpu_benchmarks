/*
Copied from: https://github.com/breagen/MachSuite/blob/master/stencil/stencil2d/stencil.c
*/

#include <stdio.h>
#include <stdlib.h>

//Define input sizes
#define col_size 64
#define row_size 128
#define f_size 9

//Data Bounds
#define TYPE int32_t
#define MAX 1000
#define MIN 1

//Set number of iterations to execute
#define MAX_ITERATION 1

void stencil( TYPE orig[row_size * col_size],
        TYPE sol[row_size * col_size],
        TYPE filter[f_size] );

struct bench_args_t {
    TYPE orig[row_size*col_size];
    TYPE sol[row_size*col_size];
    TYPE filter[f_size];
};

void stencil (TYPE orig[row_size * col_size], TYPE sol[row_size * col_size], TYPE filter[f_size]){
    int r, c, k1, k2;
    TYPE temp, mul;

    stencil_label1:for (r=0; r<row_size-2; r++) {
        stencil_label2:for (c=0; c<col_size-2; c++) {
            temp = (TYPE)0;
            stencil_label3:for (k1=0;k1<3;k1++){
                stencil_label4:for (k2=0;k2<3;k2++){
                    mul = filter[k1*3 + k2] * orig[(r+k1)*col_size + c+k2];
                    temp += mul;
                }
            }
            sol[(r*col_size) + c] = temp;
        }
    }
}
