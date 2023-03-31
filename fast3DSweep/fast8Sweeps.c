#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

double** new2DMatrix(int rows, int cols){
    
    int i, j;
    double** rval = (double**) malloc(rows*sizeof(double*));
    
    /* create 2D empty array (the idea is first create a one dimensional array of
     * pointers, and then, for each array entry, create another 1D array of pointers)*/
    for (i = 0; i < rows; i++)
        rval[i] = (double*) malloc(cols*sizeof(double));
    
    // assign values to the 2D array
    for (j = 0; j < cols; j++)
        for (i = 0; i < rows; i++)
            rval[i][j] = 0;
    return rval;
}

double*** new3DMatrix(int rows, int cols, int heis){
    
    int i,j,k;
    double*** rval = (double***) malloc(rows*sizeof(double**));
    
    /* create 3D array analogue to 2D array */
    for (i = 0; i < rows; i++){
        rval[i] = (double**) malloc(cols*sizeof(double*));
        for (j = 0; j < cols; j++){
            rval[i][j] = (double*) malloc(heis*sizeof(double));
        }
    }
    
    // assgin values to the 3D array
    for (k = 0; k < heis; k++)
        for (j = 0; j < cols; j++)
            for (i = 0; i < rows; i++)
                rval[i][j][k] = 0;
    return rval;
}

double*** get3dArray(const mxArray *mx, int rows, int cols, int heis){
    
    double* imgInPr = mxGetPr(mx);
    int i, j, k;
    double*** rval = (double***) malloc(rows*sizeof(double**));
    
    // create 3D array
    for (i = 0; i < rows; i++){
        rval[i] = (double**) malloc(cols*sizeof(double*));
        for (j = 0; j < cols; j++){
            rval[i][j] = (double*) malloc(heis*sizeof(double));
        }
    }
    
    /* assgin values to the 3D array. As array elements are column-major format
     * rval[i,j,k] = imgInPr[i+rows*(j+cols*k)] */
    for (k = 0; k < heis; k++)
        for (j = 0; j < cols; j++)
            for (i = 0; i < rows; i++)
                rval[i][j][k] = imgInPr[i+rows*(j+cols*k)];
    return rval;
}

void deleteMatrix(double ** a){
    free(a[0]);
    free(a);
}

double square(double x){
    return x * x;
}

double fastSweepingMethod(double f, double temp_b, double temp_a, double temp_c, double oldU){
    
    double a1, a2, a3, swap, x, newU;
    double *temp_all;
    int i, j, position;
    
    temp_all = (double*) malloc(3*sizeof(double));
    temp_all[0] = temp_a;
    temp_all[1] = temp_b;
    temp_all[2] = temp_c;
    
    // sort temp_all in acsending order
    for (i = 0 ; i < 2 ; i++)
    {
        position = i;
        
        for (j = i+1 ; j < 3 ; j++)
        {
            if (temp_all[position] > temp_all[j])
                position = j;
        }
        if ( position != i )
        {
            swap = temp_all[i];
            temp_all[i] = temp_all[position];
            temp_all[position] = swap;
        }
    }
    
    a1 = temp_all[0];
    a2 = temp_all[1];
    a3 = temp_all[2];
    
    // method 1 to update u
    if (square(f) >= square(a1 - a3) + square(a2 - a3)){
        x = (a1 + a2 + a3 + sqrt(3 * square(f) - square(a1 - a2) - square(a1 - a3) - square(a2 - a3))) / 3;
        newU = min(x, oldU);
    }
    else if (square(f) >= square(a1 - a2)){
        x = 0.5*(a1 + a2 + sqrt(2 * square(f) - square(a1 - a2)));
        newU = min(x, oldU);
    }
    else
        newU = min(a1 + f, oldU);
    
    return newU;
    
}

///////////////////////////////////////////////////////////////////////////
void firstSweeping(double ***weight, double ***pos, double ***u, int m, int n, int h){
    
    int i, j, k;
    double A, B, C;
    
    i = 0; // plane i = 1
    for (k = 0; k < h; k++){
        for (j = 0; j < n; j++){
            
            if (k == 0){
                if (j == 0){// at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)) { // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1) {// at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = 1; i < m - 1; i++){ // central voxels
        for (k = 0; k < h; k++){
            for (j = 0; j < n; j++){
                
                if (k == 0){
                    if (j == 0) { // at (1,1,1)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,1)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if ((k > 0) && (k < h - 1)){
                    if (j == 0) { // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on face i=1 and central pixels
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if (k == h - 1){
                    if (j == 0){// at (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
            }
        }
    }
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    i = m - 1; // plane i = 1
    for (k = 0; k < h; k++){
        for (j = 0; j < n; j++){
            
            if (k == 0){
                if (j == 0){ // at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void secondSweeping(double ***weight, double ***pos, double ***u, int m, int n, int h){
    
    int i, j, k;
    double A, B, C;
    
    ///////////////////////////////////////////////////////////////////////////
    i = m - 1; // plane i = 1
    for (k = h - 1; k >= 0; k--){
        for (j = n - 1; j >= 0; j--){
            
            if (k == 0){
                if (j == 0){ // at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = m - 2; i > 0; i--){ // central voxels
        for (k = h - 1; k >= 0; k--){
            for (j = n - 1; j >= 0; j--){
                
                if (k == 0){
                    if (j == 0) { // at (1,1,1)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,1)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if ((k > 0) && (k < h - 1)){
                    if (j == 0) { // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on face i=1 and central pixels
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if (k == h - 1){
                    if (j == 0){// at (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    i = 0; // plane i = 1
    for (k = h - 1; k >= 0; k--){
        for (j = n - 1; j >= 0; j--){
            
            if (k == 0){
                if (j == 0){// at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)) { // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1) {// at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void thirdSweeping(double ***weight, double ***pos, double ***u, int m, int n, int h){
    
    int i, j, k;
    double A, B, C;
    
    ///////////////////////////////////////////////////////////////////////////
    i = m - 1; // plane i = 1
    for (k = 0; k < h; k++){
        for (j = 0; j < n; j++){
            
            if (k == 0){
                if (j == 0){ // at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = m - 2; i > 0; i--){ // central voxels
        for (k = 0; k < h; k++){
            for (j = 0; j < n; j++){
                
                if (k == 0){
                    if (j == 0) { // at (1,1,1)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,1)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if ((k > 0) && (k < h - 1)){
                    if (j == 0) { // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on face i=1 and central pixels
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if (k == h - 1){
                    if (j == 0){// at (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    i = 0; // plane i = 1
    for (k = 0; k < h; k++){
        for (j = 0; j < n; j++){
            
            if (k == 0){
                if (j == 0){// at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)) { // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1) {// at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void fourthSweeping(double ***weight, double ***pos, double ***u, int m, int n, int h){
    
    int i, j, k;
    double A, B, C;
    
    ///////////////////////////////////////////////////////////////////////////
    i = 0; // plane i = 1
    for (k = h - 1; k >= 0; k--){
        for (j = n - 1; j >= 0; j--){
            
            if (k == 0){
                if (j == 0){// at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)) { // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1) {// at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
    
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = 1; i < m - 1; i++){ // central voxels
        for (k = h - 1; k >= 0; k--){
            for (j = n - 1; j >= 0; j--){
                
                if (k == 0){
                    if (j == 0) { // at (1,1,1)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,1)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if ((k > 0) && (k < h - 1)){
                    if (j == 0) { // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on face i=1 and central pixels
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if (k == h - 1){
                    if (j == 0){// at (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
            }
        }
    }
    
    
///////////////////////////////////////////////////////////////////////////
    i = m - 1; // plane i = 1
    for (k = h - 1; k >= 0; k--){
        for (j = n - 1; j >= 0; j--){
            
            if (k == 0){
                if (j == 0){ // at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void fifthSweeping(double ***weight, double ***pos, double ***u, int m, int n, int h){
    
    int i, j, k;
    double A, B, C;
    
    ///////////////////////////////////////////////////////////////////////////
    i = m - 1; // plane i = 1
    for (k = 0; k < h; k++){
        for (j = n - 1; j >= 0; j--){
            
            if (k == 0){
                if (j == 0){ // at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = m - 2; i > 0; i--){ // central voxels
        for (k = 0; k < h; k++){
            for (j = n - 1; j >= 0; j--){
                
                if (k == 0){
                    if (j == 0) { // at (1,1,1)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,1)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if ((k > 0) && (k < h - 1)){
                    if (j == 0) { // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on face i=1 and central pixels
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if (k == h - 1){
                    if (j == 0){// at (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    i = 0; // plane i = 1
    for (k = 0; k < h; k++){
        for (j = n - 1; j >= 0; j--){
            
            if (k == 0){
                if (j == 0){// at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)) { // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1) {// at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void sixthSweeping(double ***weight, double ***pos, double ***u, int m, int n, int h){
    
    int i, j, k;
    double A, B, C;
    
///////////////////////////////////////////////////////////////////////////
    i = 0; // plane i = 1
    for (k = h - 1; k >= 0; k--){
        for (j = 0; j < n; j++){
            
            if (k == 0){
                if (j == 0){// at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)) { // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1) {// at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
    
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = 1; i < m - 1; i++){ // central voxels
        for (k = h - 1; k >= 0; k--){
            for (j = 0; j < n; j++){
                
                if (k == 0){
                    if (j == 0) { // at (1,1,1)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,1)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if ((k > 0) && (k < h - 1)){
                    if (j == 0) { // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on face i=1 and central pixels
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if (k == h - 1){
                    if (j == 0){// at (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
            }
        }
    }
    
///////////////////////////////////////////////////////////////////////////
    i = m - 1; // plane i = 1
    for (k = h - 1; k >= 0; k--){
        for (j = 0; j < n; j++){
            
            if (k == 0){
                if (j == 0){ // at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void seventhSweeping(double ***weight, double ***pos, double ***u, int m, int n, int h){
    
    int i, j, k;
    double A, B, C;
    
    i = 0; // plane i = 1
    for (k = 0; k < h; k++){
        for (j = n - 1; j >= 0; j--){
            
            if (k == 0){
                if (j == 0){// at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)) { // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1) {// at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
    
    
    ///////////////////////////////////////////////////////////////////////////
    for (i = 1; i < m - 1; i++){ // central voxels
        for (k = 0; k < h; k++){
            for (j = n - 1; j >= 0; j--){
                
                if (k == 0){
                    if (j == 0) { // at (1,1,1)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,1)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if ((k > 0) && (k < h - 1)){
                    if (j == 0) { // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on face i=1 and central pixels
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if (k == h - 1){
                    if (j == 0){// at (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    i = m - 1;
    for (k = 0; k < h; k++){
        for (j = n - 1; j >= 0; j--){
            
            if (k == 0){
                if (j == 0){ // at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void eighthSweeping(double ***weight, double ***pos, double ***u, int m, int n, int h){
    
    int i, j, k;
    double A, B, C;
    
    ///////////////////////////////////////////////////////////////////////
    i = m - 1;
    for (k = h - 1; k >= 0; k--){
        for (j = 0; j < n; j++){
            
            if (k == 0){
                if (j == 0){ // at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i - 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////
    for (i = m - 2; i >= 1; i--){ // central voxels
        for (k = h - 1; k >= 0; k--){
            for (j = 0; j < n; j++){
                
                if (k == 0){
                    if (j == 0) { // at (1,1,1)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,1)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k + 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if ((k > 0) && (k < h - 1)){
                    if (j == 0) { // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on face i=1 and central pixels
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = min(u[i][j][k - 1], u[i][j][k + 1]);
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
                
                if (k == h - 1){
                    if (j == 0){// at (1,1,h)
                        A = u[i][j + 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if ((j > 0) && (j < n - 1)){// on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                        A = min(u[i][j - 1][k], u[i][j + 1][k]);
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                    
                    if (j == n - 1){ // at (1,n,h)
                        A = u[i][j - 1][k];
                        B = min(u[i - 1][j][k], u[i + 1][j][k]);
                        C = u[i][j][k - 1];
                        u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                        
                        if (pos[i][j][k] == 0)
                            u[i][j][k] = 0;
                    }
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////
    i = 0; // plane i = 1
    for (k = h - 1; k >= 0; k--){
        for (j = 0; j < n; j++){
            
            if (k == 0){
                if (j == 0){// at (1,1,1)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)) { // on boundary i=1 and k=1 exclude (1,1,1) and (1,n,1)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1) {// at (1,n,1)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k + 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if ((k > 0) && (k < h - 1)){
                if (j == 0){ // on boundary i=1 and j=1 excluding (1,1,1) and (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on face i=1 and central pixels
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // on boundary i=1 and j=n excluding (1,n,1) and (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = min(u[i][j][k - 1], u[i][j][k + 1]);
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
            
            if (k == h - 1){
                if (j == 0){ // at (1,1,h)
                    A = u[i][j + 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if ((j > 0) && (j < n - 1)){ // on boundary i=1 and k=h exclude (1,1,h) and (1,n,h)
                    A = min(u[i][j - 1][k], u[i][j + 1][k]);
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
                
                if (j == n - 1){ // at (1,n,h)
                    A = u[i][j - 1][k];
                    B = u[i + 1][j][k];
                    C = u[i][j][k - 1];
                    u[i][j][k] = fastSweepingMethod(weight[i][j][k], A, B, C, u[i][j][k]);
                    
                    if (pos[i][j][k] == 0)
                        u[i][j][k] = 0;
                }
            }
        }
    }
}

// main function using matlab mex interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    int rows, cols, heis, i, j, k;
    double ***weight, ***pos, ***u, *imgOutPr;
    
    /* Size of input */
    int nsubs = 0;
    const mwSize *idims;
    
    /*  Get the number of dimensions */
    nsubs = mxGetNumberOfDimensions(prhs[0]);
    
    /* Get the sizes of the inputs */
    idims = mxGetDimensions(prhs[0]);
    rows = idims[0];
    cols = idims[1];
    heis = idims[2];
    
    /* transform the input 1D array to 3D array */
    weight = get3dArray(prhs[0], rows, cols, heis);
    pos = get3dArray(prhs[1], rows, cols, heis);
    u = get3dArray(prhs[2], rows, cols, heis);
    
    firstSweeping(weight, pos, u, rows, cols, heis);
    secondSweeping(weight, pos, u, rows, cols, heis);
    thirdSweeping(weight, pos, u, rows, cols, heis);
    fourthSweeping(weight, pos, u, rows, cols, heis);
    fifthSweeping(weight, pos, u, rows, cols, heis);
    sixthSweeping(weight, pos, u, rows, cols, heis);
    seventhSweeping(weight, pos, u, rows, cols, heis);
    eighthSweeping(weight, pos, u, rows, cols, heis);
    
    // output the the matrix
    plhs[0] = mxCreateNumericArray(nsubs, idims, mxDOUBLE_CLASS, mxREAL);
    imgOutPr = mxGetPr(plhs[0]);
    
    for (k = 0; k < heis; k++)
        for (j = 0; j < cols; j++)
            for (i = 0; i < rows; i++)
                imgOutPr[i + rows*(j + cols*k)] = u[i][j][k];
    
    return;
}