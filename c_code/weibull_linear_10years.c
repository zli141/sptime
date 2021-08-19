
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "R_ext/Applic.h"
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>        /* for Lapack (dpotrf, etc.) and BLAS */
#include <R_ext/Linpack.h>


double **create_matrixd (int row, int col)
{
    register int i;
    double **matrix;
    
    matrix = Calloc(row, double *);
    assert(matrix);
    for(i = 0; i < row; ++i) {
        matrix[i] = Calloc(col, double);
        assert(matrix[i]);
    }
    return matrix;
}

double *create_vectord (int dim)
{
    double *vector;
    
    vector = Calloc (dim, double);
    assert(vector);
    
    return vector;
}


int *create_vector (int dim)
{
    int *vector;
    
    vector = Calloc (dim, int);
    assert(vector);
    
    return vector;
}


int **create_matrix (int row, int col)
{
    register int i;
    int **matrix;
    
    matrix = Calloc(row, int *);
    assert(matrix);
    for(i = 0; i < row; ++i) {
        matrix[i] = Calloc(col,int);
        assert(matrix[i]);
    }
    return matrix;
    
}

void matvector (double *vector, double **matrix,int nrow, int ncol)
{
    int i, j;
    
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) vector[j*nrow+i]=matrix[i][j];
    }
    
}

void matsclpr (double **matrix, double a, double **b, int rows, int cols)
{
    int i, j;
    
   	for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            matrix[i][j] = a * b[i][j];
        }
   	}
}

void matsum (double **matrix, double **a, double **b, int rows, int cols)
{
    int i, j;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) matrix[i][j] = a[i][j] + b[i][j];
    }
}


void matdiff (double **matrix, double **a, double **b, int rows, int cols)
{
    int i, j;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) matrix[i][j] = a[i][j] - b[i][j];
    }
}


void freematd (double **a, int row, int col)
{
    int i;
    
    for (i = 0; i < row; i++) Free(a[i]);
    Free(a);
}

void zero_mat(double **matrix, int rows, int cols)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++) 	matrix[i][j]=0;
    }
    
}
void ident_mat(double **matrix, int dim)
{
    int i, j;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        { if (i!=j) matrix[i][j]=0;
	       else matrix[i][j]=1;
        }
    }
    
}

void vectormat (double **matrix, double *vector, int nrow, int ncol)
{
    int i, j;
    
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) matrix[i][j] = vector[j*nrow+i];
    }
    
}


double *inverse (double *invertmat, int dim, double **matrix)
{
    
    double *qraux, *b,  *y, **ymat, *temp1, *temp2;
    double *tol, tol1;
    int *pivot, *n, *k, rank, *inf, l;
    n=&dim;
    pivot=create_vector(dim);
    qraux=create_vectord(dim);
    b=create_vectord(dim*2);
    k=&rank;
    tol1=1.0e-07;
    tol=&tol1;
    inf=&l;
    
    temp1=create_vectord(dim*dim);
    matvector (temp1, matrix,dim, dim);
    
    ymat=create_matrixd(dim, dim);
    ident_mat(ymat, dim);
    y=create_vectord(dim*dim);
    matvector (y, ymat,dim, dim);
    
    temp2=create_vectord(dim*dim);
    
    F77_CALL(dqrdc2)(temp1,n,n,n,tol,k,qraux,pivot,b);
    
   
    
    F77_CALL(dqrcf)(temp1,n,k,qraux,y, n, temp2,inf);
   
    freematd(ymat,dim,dim);
    Free(pivot);
    Free(qraux);
    Free(b);
    Free(temp1);
    Free(y);   
    return temp2;
    
}


void transpose (double **matrix, double **x, int rows, int cols)
{
    int  i, j;
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++) {
            matrix[j][i] = x[i][j];
        }
    }
}

void matvecpr (double *vector, double **a, int rows, int cols, double *x)
{
    int i, j;
    double sum;
    
    for (i = 0; i < rows; i++) {
        sum = 0;
        for (j = 0; j < cols; j++) {
            sum += a[i][j] * x[j];
        }
        vector[i] = sum;
    }
}

void matpr (double **c, double **a, double **b, int m, int p, int n)
{
    int i, j, k;
    double sum;
    
    for (i = 0; i < m ; i++) {
        for(j = 0; j < n; j++) {
            sum = 0 ;
            for (k = 0; k < p; k++) {
                sum += a[i][k] * b[k][j];
            }
            c[i][j]  = sum;
        }
    }
}


double *std_rWishart_factor(double nu, int p, int upper, double ans[])
{
    int pp1 = p + 1;
    memset(ans, 0, p * p * sizeof(double));
    
    for (int j = 0; j < p; j++) {	/* jth column */
        ans[j * pp1] = sqrt(rchisq(nu - (double) j));
        for (int i = 0; i < j; i++) {
            int uind = i + j * p, /* upper triangle index */
            lind = j + i * p; /* lower triangle index */
            ans[(upper ? uind : lind)] = norm_rand();
            ans[(upper ? lind : uind)] = 0;
        }
    }
    return ans;
}

void rwishart(double *ansj, double nuP, double *scCp, int *dims){
    int info,  psqr;
    double *tmp, nu = nuP, one = 1, zero = 0;
    // allocate early to avoid memory leaks in Callocs below.
    psqr = dims[0] * dims[0];
    tmp = Calloc(psqr, double);
    F77_CALL(dpotrf)("U", &(dims[0]), scCp, &(dims[0]), &info);
    
    std_rWishart_factor(nu, dims[0], 1, tmp);
    F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims,
                    &one, scCp, dims, tmp, dims);
    F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
                    &one, tmp, &(dims[1]),
                    &zero, ansj, &(dims[1]));
    
    for (int i = 1; i < dims[0]; i++){
        for (int k = 0; k < i; k++){
            ansj[i + k * dims[0]] = ansj[k + i * dims[0]];
        }
    }
    Free(tmp);
}


/*********beta0 *****************/

double target_b (double b0, double b1, double b2, double b3,double b14, double b15, double b16, double b17, double b18, double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2, double sig, int *Cores, double *year, double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
    for (int i = 0; i < N; i++)
    {
        mu=b0 + b1*x1[i]+b2*x2[i]+b3*(year[i]-1)+b14*x3[i]+b15*x4[i]+b16*x5[i]+b17*x6[i]+b18*x7[i]+b19*x8[i]+b20*x9[i]+b21*x10[i]+b22*x11[i]+b23*x12[i]+phi1[Cores[i]]+phi2[Cores[i]]*(year[i]-1);
        lam=exp(-mu/sig);
        a=1/(sig);
        sum=sum+s[i]*log(lam*a*pow(y[i],(a-1)))-lam*pow(y[i],a);
        
    }
        return sum;
}


/*********for phi *****************/

double target_phi(double b0, double b1, double b2, double b3, double b14,double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2, double sig,int *Cores, double *year, double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N, int no_regions, double **RTR11P1, double **lamda, double x1_new, double x2_new, int k)
{
    double sum = 0.0, func=0.0,mu=0.0,a=0.0,lam=0.0;
    int i ,j, m;
    
    double *temp1, *temp2;
    temp1=create_vectord(no_regions);
    temp2=create_vectord(no_regions);
    
    
    for(m=0; m<no_regions; m++)
    {
      if(m==k) {
        temp1[m]=x1_new;
        temp2[m]=x2_new;
       }
    else
    {
        temp1[m]=phi1[m];
        temp2[m]=phi2[m];
     }
    }

    
    for (i = 0; i < N; i++)
    {
        mu=b0+b1*x1[i]+b2*x2[i]+b14*x3[i]+b15*x4[i]+b16*x5[i]+b17*x6[i]+b18*x7[i]+b19*x8[i]+b20*x9[i]+b21*x10[i]+b22*x11[i]+b23*x12[i]+b3*(year[i]-1)+temp1[Cores[i]]+temp2[Cores[i]]*(year[i]-1);
        lam=exp(-mu/sig);
        a=1/(sig);
        sum=sum+s[i]*log(lam*a*pow(y[i],(a-1)))-lam*pow(y[i],a);
        
    }
    
    
    func=func+lamda[0][0]* RTR11P1[k][k]*temp1[k]*temp1[k]+lamda[1][1]* RTR11P1[k][k]*temp2[k]*temp2[k]+2*RTR11P1[k][k]*temp1[k]*temp2[k]*lamda[0][1];
    
    for (j=0; j<no_regions; j++)
    { if (j != k)
		  {
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp1[k]*lamda[0][0]+2*RTR11P1[j][k]*temp2[j]*temp2[k]*lamda[1][1];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp2[k]*lamda[0][1]+2*RTR11P1[k][j]*temp1[k]*temp2[j]*lamda[0][1];
          }
    }
    sum=sum-0.5*func;
    Free(temp1);
    Free(temp2);
    return (sum);
}

double target_sig (double b0, double b1, double b2, double b3,double b14, double b15, double b16, double b17, double b18, double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2, double x, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
     a=1/(x);
  
    for (int i = 0; i < N; i++)
    {
        mu=b0 + b1*x1[i]+b2*x2[i]+b3*(year[i]-1)+b14*x3[i]+b15*x4[i]+b16*x5[i]+b17*x6[i]+b18*x7[i]+b19*x8[i]+b20*x9[i]+b21*x10[i]+b22*x11[i]+b23*x12[i]+phi1[Cores[i]]+phi2[Cores[i]]*(year[i]-1);
        lam=exp(-mu/x);
        sum=sum+s[i]*log(lam*a*pow(y[i],(a-1)))-lam*pow(y[i],a);
      
        
    }
    
    sum=sum+dgamma(x,0.001,1000,1);
    
    return sum;
}


/***** for lamda matrix ***********************/
void target_lamda(double **matrix, double *phi1, double *phi2, int no_regions, double **RTR11P1, int p)
{
    double **temp;
    int i, j;
    temp=create_matrixd(p,p);
    zero_mat(temp, p, p);
    
    
    for( i=0; i<no_regions; i++)
    {
        for (j=0; j<no_regions; j++)
    
        {
        temp[0][0]=temp[0][0]+ phi1[i]*phi1[j]*RTR11P1[i][j];
        temp[0][1]=temp[0][1]+ phi1[i]*phi2[j]*RTR11P1[i][j];
        temp[1][1]=temp[1][1]+ phi2[i]*phi2[j]*RTR11P1[i][j];
        }
        
    }
    temp[1][0]=temp[0][1];
    transpose(matrix, temp, p, p);
    
    freematd(temp,p, p);
    
}


/********* Metroplis-Hasting for rho and beta0  ***********************/

double metroplis_b0(double x1_prev, double b1,double b2,double b3,double b14, double b15, double b16, double b17, double b18, double b19,  double b20, double b21, double b22, double b23, double *phi1, double *phi2, double sig,  int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);

    y1=x1_prev+0.1*rnorm(0,1);
  
    temp1=target_b(y1,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(x1_prev,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=x1_prev;}
    }
    else {newx=x1_prev;}

    return (newx);
}


double metroplis_b1(double b0, double x_prev,double b2,double b3, double b14, double b15, double b16, double b17, double b18, double b19,  double b20, double b21, double b22, double b23, double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
  
    temp1=target_b(b0,y_new,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,x_prev,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    
    printf("v1 y_new is %f \n",y_new);   
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}


double metroplis_b2( double b0, double b1, double x_prev,double b3, double b14, double b15, double b16, double b17, double b18, double b19, double b20, double b21, double b22, double b23, double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,y_new,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,x_prev,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b3( double b0, double b1,double b2,double x_prev,double b14, double b15, double b16, double b17, double b18, double b19,double b20, double b21, double b22, double b23,  double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,y_new,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,x_prev,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b14( double b0, double b1,double b2,double b3,double x_prev,double b15, double b16, double b17, double b18, double b19,double b20, double b21, double b22, double b23, double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,y_new,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,x_prev,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b15( double b0, double b1,double b2,double b3, double b14,double x_prev, double b16, double b17, double b18, double b19,double b20, double b21, double b22, double b23, double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,y_new,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,x_prev,b16,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b16( double b0, double b1,double b2,double b3, double b14, double b15,double x_prev, double b17, double b18, double b19, double b20, double b21, double b22, double b23, double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,b15,y_new,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,b15,x_prev,b17,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b17( double b0, double b1,double b2,double b3, double b14, double b15,double b16, double x_prev, double b18, double b19,double b20, double b21, double b22, double b23, double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,b15,b16,y_new,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,b15,b16,x_prev,b18,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b18( double b0, double b1,double b2,double b3, double b14, double b15,double b16, double b17,double x_prev, double b19,double b20, double b21, double b22, double b23, double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,b15,b16,b17,y_new,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,b15,b16,b17,x_prev,b19,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b19( double b0, double b1,double b2,double b3, double b14, double b15,double b16, double b17, double b18, double x_prev, double b20, double b21, double b22, double b23,double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,y_new,b20,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,x_prev,b20,b21,b22,b23, phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b20( double b0, double b1,double b2,double b3, double b14, double b15,double b16, double b17, double b18, double b19, double x_prev, double b21, double b22, double b23,double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,y_new,b21,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,x_prev,b21,b22,b23, phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b21( double b0, double b1,double b2,double b3, double b14, double b15,double b16, double b17, double b18, double b19, double b20, double x_prev, double b22, double b23,double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,y_new,b22,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,x_prev,b22,b23, phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b22( double b0, double b1,double b2,double b3, double b14, double b15,double b16, double b17, double b18, double b19, double b20, double b21, double x_prev, double b23, double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,y_new,b23,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,x_prev,b23, phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double metroplis_b23( double b0, double b1,double b2,double b3, double b14, double b15,double b16, double b17, double b18, double b19, double b20, double b21, double b22, double x_prev,double *phi1,double *phi2, double sig, int *Cores, double *year,double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,y_new,phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,x_prev, phi1,phi2,sig,Cores, year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}

double *center(double *x,int n_region){
    double sum=0.0, avg;
    for (int i=0;i<n_region;i++){
        sum=sum+x[i];
    }
    avg=sum/n_region;
    for (int j=0;j<n_region;j++){
        x[j]=x[j]-avg;
    }
    return (x);
}

/********* Metroplis-Hasting for random effects phi ***********************/

void metroplis_phi(double b0, double b1,double b2,double b3,double b14, double b15, double b16, double b17, double b18, double b19, double b20, double b21, double b22, double b23, double *phi_curr1, double *phi_curr2, double sig, double **lamda, int *Cores, double *year, double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y,int *s, int N, int no_regions, double **RTR11P1,double x1_prev, double x2_prev, int j)
{
    double lgratio=0, u ,temp1=0, temp2=0;
    double x1_curr, x2_curr;
   
    
    u=runif(0,1);
    
    x1_curr=x1_prev+rnorm(0,0.1);
    x2_curr=x2_prev+rnorm(0,0.1);
    
    temp1=target_phi(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi_curr1,phi_curr2,sig, Cores, year, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N, no_regions,RTR11P1,lamda, x1_curr, x2_curr,j);
    temp2=target_phi(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi_curr1,phi_curr2,sig, Cores, year, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N, no_regions,RTR11P1,lamda, x1_prev, x2_prev,j);
    lgratio=temp1-temp2;
  
    if (u>0)
    { if (lgratio>=log(u))
    {   phi_curr1[j]=x1_curr;
        phi_curr2[j]=x2_curr;
        
    }
    else
    {
        phi_curr1[j]=x1_prev;
        phi_curr2[j]=x2_prev;
        
    }
    }
    else
        
    {   phi_curr1[j]=x1_prev;
        phi_curr2[j]=x2_prev;
        
    }
    
}



double metroplis_sig( double b0, double b1, double b2, double b3, double b14, double b15, double b16, double b17, double b18, double b19, double b20, double b21, double b22, double b23,double *phi_curr1, double *phi_curr2, double x_prev, int *Cores, double *year,double *x1, double *x2, double *x3, double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12, double *y,int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+rnorm(0,0.1);
    if (y_new>0){
    temp1=target_sig(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi_curr1, phi_curr2,y_new,Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_sig(b0,b1,b2,b3,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi_curr1,phi_curr2,x_prev,Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    }
    else x=x_prev;
    
    
    
    return(x);
}


double LL (double b0, double b1, double b2, double b3, double b14, double b15, double b16, double b17, double b18, double b19,double b20, 
              double b21, double b22, double b23, double *phi1, double *phi2, double sig, int *Cores, double *year, double *x1, double *x2,double *x3, 
                 double *x4,double *x5, double *x6,double *x7, double *x8,double *x9, double *x10, double *x11, double *x12,  double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
    for (int i = 0; i < N; i++){
        mu=b0 + b1*x1[i]+b2*x2[i]+b3*(year[i]-1)+b14*x3[i]+b15*x4[i]+b16*x5[i]+b17*x6[i]+b18*x7[i]+b19*x8[i]+b20*x9[i]+b21*x10[i]+b22*x11[i]+b23*x12[i]+phi1[Cores[i]]+phi2[Cores[i]]*(year[i]-1);
        lam=exp(-mu/sig);
        a=1/(sig);
        sum=sum+s[i]*log(lam*a*pow(y[i],(a-1)))-lam*pow(y[i],a);
    }
    return -2*sum;
}


void weibul(double *b0, double *b1, double *b2, double *b3,double *b14, double *b15, double *b16, double *b17, double *b18, double *b19, double *b20, 
           double *b21, double *b22, double *b23,  double *sig, double *lamda_p0,
           double *phi1, double *phi2, double *vecD, double *vecC,double *y, int *s, int *N_p, 
       double *x1, double *x2,double *x3, double *x4,double *x5, double *x6,double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,  
       double *year, int *Cores, int *nend, int *no_regions_p,int *dims,double *lamda_prec,double *sumLnF,double *phi_sample,int *n_burn){
	int i, j;
    int p=dims[0];
    int N=N_p[0];
    int no_regions=no_regions_p[0];
    double **temp_lamda, **tempsum, *tempsum_INV, **lamda, **lamda_prec_mat,**D,**C,**RTR11P1;
    double x1_prev,x2_prev;
	temp_lamda=create_matrixd(p,p);
    tempsum=create_matrixd(p,p);
	tempsum_INV=create_vectord(p*p);
    lamda=create_matrixd(p,p);
	lamda_prec_mat=create_matrixd(p,p);
    D=create_matrixd(no_regions,no_regions);
    C=create_matrixd(no_regions,no_regions);
    RTR11P1=create_matrixd(no_regions,no_regions);
	
    vectormat(lamda, lamda_p0, p, p);//transfer lamda_p0 the init into init matrix lamda
    vectormat(lamda_prec_mat,lamda_prec,p,p);//transfer wishart init precision matrix
    
    vectormat(D,vecD,no_regions,no_regions);//transfer diagnol matrix
    vectormat(C,vecC,no_regions,no_regions);//transfer adjacent matrix
    
    matdiff(RTR11P1, D, C, no_regions, no_regions);//D_w-C

	 GetRNGstate( );
    
   for (i = 1; i < nend[0]; i++)
     //for (i = 1; i < 2; i++)
	{
    
	 b0[i]=metroplis_b0(b0[i-1],b1[i-1],b2[i-1],b3[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
     b1[i]=metroplis_b1(b0[i],b1[i-1],b2[i-1],b3[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
     b2[i]=metroplis_b2(b0[i],b1[i],b2[i-1],b3[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
     b3[i]=metroplis_b3(b0[i],b1[i],b2[i],b3[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
     b14[i]=metroplis_b14(b0[i],b1[i],b2[i],b3[i],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
     b15[i]=metroplis_b15(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b16[i]=metroplis_b16(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b17[i]=metroplis_b17(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i-1],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b18[i]=metroplis_b18(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i-1],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b19[i]=metroplis_b19(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i-1],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b20[i]=metroplis_b20(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i-1], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
     b21[i]=metroplis_b21(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i], b21[i-1], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
   b22[i]=metroplis_b22(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i], b21[i], b22[i-1],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b23[i]=metroplis_b23(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i], b21[i], b22[i],b23[i-1],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    
    
    
    for (j = 0; j < no_regions; j++)
	    {
            x1_prev=phi1[j];
            x2_prev=phi2[j];
            metroplis_phi(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i], b21[i], b22[i],b23[i],phi1,phi2,sig[i-1],lamda, Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N,no_regions,RTR11P1,x1_prev,x2_prev,j);
        }
    //printf("Before center %d is %f %f %f \n",i+1,phi1[0],phi1[1],phi1[2]);
    //printf("Before center %d is %f %f %f \n",i+1,phi2[0],phi2[1],phi2[2]);
      phi1=center(phi1,no_regions);
      phi2=center(phi2,no_regions);
      sig[i]=metroplis_sig(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i], b21[i], b22[i],b23[i],phi1,phi2,sig[i-1],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);

      target_lamda(temp_lamda,phi1,phi2,no_regions,RTR11P1,p);
      matsum (tempsum, lamda_prec_mat, temp_lamda, p, p);
      tempsum_INV=inverse(tempsum_INV,p,tempsum);
      rwishart(lamda_p0,p+no_regions,tempsum_INV,dims);
      vectormat(lamda,lamda_p0,p,p);
      if(i>=n_burn[0]) sumLnF[0]=sumLnF[0]+LL(b0[i],b1[i],b2[i],b3[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i], b21[i], b22[i],b23[i-1],phi1,phi2,sig[i],Cores,year,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        for(int k=0;k<no_regions;k++){
          phi_sample[k*nend[0]+i]=phi1[k];
          phi_sample[(k+no_regions)*nend[0]+i]=phi2[k];
        }
    }
    
    freematd(tempsum,p, p);
    Free(tempsum_INV);
    freematd(temp_lamda, p, p);
    freematd(lamda_prec_mat, p, p);
    freematd(RTR11P1,no_regions,no_regions);
    freematd(C, no_regions, no_regions);
    freematd(D, no_regions, no_regions);
    freematd(lamda, p,p);
    PutRNGstate( );

}






