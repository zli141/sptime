
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

/*void multinormal (double *vector, double *mean, double **disp, int dim)
{
    int i, *inf, l, *p;
    double *x, *z, **A, *temp1, *temp2, **tA  ;
    
    temp1= create_vectord(dim*dim);
    temp2= create_vectord (dim*dim);
    p=&dim; inf=&l;
    matvector (temp1, disp,dim, dim);
    F77_CALL(dpotrf)("U", p, temp1, p, 0);
    //F77_CALL(chol)(temp1, p, p, temp2, inf);
    if(inf[0]!=0)
    {
        REprintf("This is not P.d matrix\n");
        Free (temp1);
        Free (temp2);
        return;
    }
    x = create_vectord (dim);
    z = create_vectord (dim);
    A=create_matrixd (dim, dim);
    tA=create_matrixd (dim, dim);
    vectormat (A, temp1, dim,dim);
    transpose (tA, A, dim, dim);
    for (i = 0; i < dim; i++) z[i] = norm_rand( );
    matvecpr (x, tA, dim, dim, z);
    for (i = 0; i < dim; i++) vector[i] = mean[i] + x[i];
    Free (x);
    Free (z);
    Free (temp1);
    Free (temp2);
    freematd (A, dim, dim);
    freematd (tA, dim, dim);
}
*/
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

double target_b (double b0, double b1, double b2, double b3, double b4, double b5, double b6, double b7, double b8, double b9, double b10, double b11,
                 double b12,double b14, double b15,double b16,double b17, double b18,double b19, double b20, double b21, double b22, 
                 double b23, double *phi1, double *phi2, double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, 
                 double *phi10,double *phi11,  double sig, int *Cores, double *year1, double *year2,double *year3,double *year4, 
                 double *year5, double *year6, double *year7, double *year8, double *year9,double *year10,  double *x1, double *x2, 
                 double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double sum=0.0,mu=0.0,lam=0.0,a= 0.0;
    for (int i = 0; i <N; i++)
    {
        mu=b0 + b1*x1[i]+b2*x2[i]+b14*x3[i]+b15*x4[i]+b16*x5[i]+b17*x6[i]+b18*x7[i]+b19*x8[i]+b20*x9[i]+b21*x10[i]+b22*x11[i]+b23*x12[i]
           +b3*year1[i]+b4*year2[i]+b5*year3[i]+b6*year4[i]+b7*year5[i]+b8*year6[i]+b9*year7[i]+b10*year8[i]+b11*year9[i]+b12*year10[i]
           +phi1[Cores[i]]+phi2[Cores[i]]*year1[i]+phi3[Cores[i]]*year2[i]+phi4[Cores[i]]*year3[i]+phi5[Cores[i]]*year4[i]+phi6[Cores[i]]*year5[i]
           +phi7[Cores[i]]*year6[i]+phi8[Cores[i]]*year7[i]+phi9[Cores[i]]*year8[i]+phi10[Cores[i]]*year9[i]+phi11[Cores[i]]*year10[i];
        lam=exp(-mu/sig);
        a=1/(sig);
        sum=sum+s[i]*log(lam*a*pow(y[i],(a-1)))-lam*pow(y[i],a);
        
    }
    
       return sum;
}





/*********for phi *****************/

double target_phi(double b0, double b1, double b2, double b3, double b4, double b5, double b6, double b7, double b8, double b9, double b10, double b11,
                  double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23,
                  double *phi1, double *phi2, double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, 
                  double *phi10,double *phi11,  double sig, int *Cores, double *year1, double *year2,double *year3,double *year4, double *year5, double *year6, 
                  double *year7, double *year8, double *year9,double *year10,  double *x1, double *x2, double  *x3, double *x4, double *x5, double *x6,
                  double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N, int no_regions, double **RTR11P1, double **lamda, 
                  double x1_new, double x2_new, double x3_new, double x4_new, double x5_new, double x6_new, double x7_new,
                   double x8_new, double x9_new, double x10_new, double x11_new, int k)
{
    double sum = 0.0, func=0.0,mu=0.0,a=0.0,lam=0.0;
    int i ,j, m;
    
    double *temp1, *temp2, *temp3, *temp4, *temp5, *temp6, *temp7, *temp8, *temp9, *temp10, *temp11;
    temp1=create_vectord(no_regions);
    temp2=create_vectord(no_regions);
    temp3=create_vectord(no_regions);
    temp4=create_vectord(no_regions);
    temp5=create_vectord(no_regions);
    temp6=create_vectord(no_regions);
   temp7=create_vectord(no_regions);
   temp8=create_vectord(no_regions);
   temp9=create_vectord(no_regions);
   temp10=create_vectord(no_regions);
   temp11=create_vectord(no_regions);
  
      
    for(m=0; m<no_regions; m++)
    {
      if(m==k) {
        temp1[m]=x1_new;
        temp2[m]=x2_new;
        temp3[m]=x3_new;
        temp4[m]=x4_new;
        temp5[m]=x5_new;
        temp6[m]=x6_new;
        temp7[m]=x7_new;
        temp8[m]=x8_new;
        temp9[m]=x9_new;
        temp10[m]=x10_new;
        temp11[m]=x11_new;  
        
    
       }
    else
    {
        temp1[m]=phi1[m];
        temp2[m]=phi2[m];
        temp3[m]=phi3[m];
        temp4[m]=phi4[m];
        temp5[m]=phi5[m];
        temp6[m]=phi6[m];
        temp7[m]=phi7[m];
        temp8[m]=phi8[m];
        temp9[m]=phi9[m];
        temp10[m]=phi10[m];
        temp11[m]=phi11[m];  
     
       
     }
    }

    
    for (i = 0; i < N; i++)
    {
         mu=b0 + b1*x1[i]+b2*x2[i]+b14*x3[i]+b15*x4[i]+b16*x5[i]+b17*x6[i]+b18*x7[i]+b19*x8[i]+b20*x9[i]+b21*x10[i]+b22*x11[i]+b23*x12[i]
           +b3*year1[i]+b4*year2[i]+b5*year3[i]+b6*year4[i]+b7*year5[i]+b8*year6[i]+b9*year7[i]+b10*year8[i]+b11*year9[i]+b12*year10[i]
           +temp1[Cores[i]]+temp2[Cores[i]]*year1[i]+temp3[Cores[i]]*year2[i]+temp4[Cores[i]]*year3[i]+temp5[Cores[i]]*year4[i]+temp6[Cores[i]]*year5[i]
           +temp7[Cores[i]]*year6[i]+temp8[Cores[i]]*year7[i]+temp9[Cores[i]]*year8[i]+temp10[Cores[i]]*year9[i]+temp11[Cores[i]]*year10[i];
        lam=exp(-mu/sig);
        a=1/(sig);
        sum=sum+s[i]*log(lam*a*pow(y[i],(a-1)))-lam*pow(y[i],a);
        
    }
    

    
  func=func
  +lamda[0][0]* RTR11P1[k][k]*temp1[k]*temp1[k]+lamda[1][1]* RTR11P1[k][k]*temp2[k]*temp2[k]
  +lamda[2][2]* RTR11P1[k][k]*temp3[k]*temp3[k]+lamda[3][3]* RTR11P1[k][k]*temp4[k]*temp4[k]
  +lamda[4][4]* RTR11P1[k][k]*temp5[k]*temp5[k]+lamda[5][5]* RTR11P1[k][k]*temp6[k]*temp6[k]
  +lamda[6][6]* RTR11P1[k][k]*temp7[k]*temp7[k]+lamda[7][7]* RTR11P1[k][k]*temp8[k]*temp8[k]
  +lamda[8][8]* RTR11P1[k][k]*temp9[k]*temp9[k]+lamda[9][9]* RTR11P1[k][k]*temp10[k]*temp10[k]
  +lamda[10][10]* RTR11P1[k][k]*temp11[k]*temp11[k]
  
  +2*RTR11P1[k][k]*temp1[k]*temp2[k]*lamda[0][1]+2*RTR11P1[k][k]*temp1[k]*temp3[k]*lamda[0][2]
  +2*RTR11P1[k][k]*temp1[k]*temp4[k]*lamda[0][3]+2*RTR11P1[k][k]*temp1[k]*temp5[k]*lamda[0][4]
  +2*RTR11P1[k][k]*temp1[k]*temp6[k]*lamda[0][5]+2*RTR11P1[k][k]*temp1[k]*temp7[k]*lamda[0][6]
  +2*RTR11P1[k][k]*temp1[k]*temp8[k]*lamda[0][7]+2*RTR11P1[k][k]*temp1[k]*temp9[k]*lamda[0][8]
  +2*RTR11P1[k][k]*temp1[k]*temp10[k]*lamda[0][9]+2*RTR11P1[k][k]*temp1[k]*temp11[k]*lamda[0][10]

  
  +2*RTR11P1[k][k]*temp2[k]*temp3[k]*lamda[1][2]+2*RTR11P1[k][k]*temp2[k]*temp4[k]*lamda[1][3]
  +2*RTR11P1[k][k]*temp2[k]*temp5[k]*lamda[1][4]+2*RTR11P1[k][k]*temp2[k]*temp6[k]*lamda[1][5]
  +2*RTR11P1[k][k]*temp2[k]*temp7[k]*lamda[1][6]+2*RTR11P1[k][k]*temp2[k]*temp8[k]*lamda[1][7]
  +2*RTR11P1[k][k]*temp2[k]*temp9[k]*lamda[1][8]+2*RTR11P1[k][k]*temp2[k]*temp10[k]*lamda[1][9]
  +2*RTR11P1[k][k]*temp2[k]*temp11[k]*lamda[1][10]
 
  +2*RTR11P1[k][k]*temp3[k]*temp4[k]*lamda[2][3]+2*RTR11P1[k][k]*temp3[k]*temp5[k]*lamda[2][4]+2*RTR11P1[k][k]*temp3[k]*temp6[k]*lamda[1][5]
  +2*RTR11P1[k][k]*temp3[k]*temp7[k]*lamda[2][6]+2*RTR11P1[k][k]*temp3[k]*temp8[k]*lamda[2][7]
  +2*RTR11P1[k][k]*temp3[k]*temp9[k]*lamda[2][8]+2*RTR11P1[k][k]*temp3[k]*temp10[k]*lamda[2][9]
  +2*RTR11P1[k][k]*temp3[k]*temp11[k]*lamda[2][10]

  +2*RTR11P1[k][k]*temp4[k]*temp5[k]*lamda[3][4]+2*RTR11P1[k][k]*temp4[k]*temp6[k]*lamda[3][5]
  +2*RTR11P1[k][k]*temp4[k]*temp7[k]*lamda[3][6]+2*RTR11P1[k][k]*temp4[k]*temp8[k]*lamda[3][7]
  +2*RTR11P1[k][k]*temp4[k]*temp9[k]*lamda[3][8]+2*RTR11P1[k][k]*temp4[k]*temp10[k]*lamda[3][9]
  +2*RTR11P1[k][k]*temp4[k]*temp11[k]*lamda[3][10]

  +2*RTR11P1[k][k]*temp5[k]*temp6[k]*lamda[4][5] +2*RTR11P1[k][k]*temp5[k]*temp7[k]*lamda[4][6]+2*RTR11P1[k][k]*temp4[k]*temp8[k]*lamda[4][7]
  +2*RTR11P1[k][k]*temp5[k]*temp9[k]*lamda[4][8]+2*RTR11P1[k][k]*temp5[k]*temp10[k]*lamda[4][9]
  +2*RTR11P1[k][k]*temp5[k]*temp11[k]*lamda[4][10]

  +2*RTR11P1[k][k]*temp6[k]*temp7[k]*lamda[5][6]+2*RTR11P1[k][k]*temp6[k]*temp8[k]*lamda[5][7]
  +2*RTR11P1[k][k]*temp6[k]*temp9[k]*lamda[5][8]+2*RTR11P1[k][k]*temp6[k]*temp10[k]*lamda[5][9]
  +2*RTR11P1[k][k]*temp6[k]*temp11[k]*lamda[5][10]

  +2*RTR11P1[k][k]*temp7[k]*temp8[k]*lamda[6][7]
  +2*RTR11P1[k][k]*temp7[k]*temp9[k]*lamda[6][8]+2*RTR11P1[k][k]*temp7[k]*temp10[k]*lamda[6][9]
  +2*RTR11P1[k][k]*temp7[k]*temp11[k]*lamda[6][10]
  
  +2*RTR11P1[k][k]*temp8[k]*temp9[k]*lamda[7][8]+2*RTR11P1[k][k]*temp8[k]*temp10[k]*lamda[7][9]
  +2*RTR11P1[k][k]*temp8[k]*temp11[k]*lamda[7][10]

  +2*RTR11P1[k][k]*temp9[k]*temp10[k]*lamda[8][9]
  +2*RTR11P1[k][k]*temp9[k]*temp11[k]*lamda[8][10]
     
  +2*RTR11P1[k][k]*temp10[k]*temp11[k]*lamda[9][10];
   

    for (j=0; j<no_regions; j++)
    { if (j != k)
		  {
             func=func+ 2*RTR11P1[j][k]*temp1[j]*temp1[k]*lamda[0][0]+2*RTR11P1[j][k]*temp2[j]*temp2[k]*lamda[1][1]
                   +2*lamda[2][2]* RTR11P1[k][k]*temp3[j]*temp3[k] +2*lamda[3][3]* RTR11P1[k][k]*temp4[j]*temp4[k]
                   +2*lamda[4][4]* RTR11P1[k][k]*temp5[j]*temp5[k] +2*lamda[5][5]* RTR11P1[k][k]*temp6[j]*temp6[k]
                   +2*lamda[6][6]* RTR11P1[k][k]*temp7[j]*temp7[k] +2*lamda[7][7]* RTR11P1[k][k]*temp8[j]*temp8[k]
                   +2*lamda[8][8]* RTR11P1[k][k]*temp9[j]*temp9[k] +2*lamda[9][9]* RTR11P1[k][k]*temp10[j]*temp10[k]
                   +2*lamda[10][10]* RTR11P1[k][k]*temp11[j]*temp11[k];
              
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp2[k]*lamda[0][1]+2*RTR11P1[k][j]*temp1[k]*temp2[j]*lamda[0][1];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp3[k]*lamda[0][2]+2*RTR11P1[k][j]*temp1[k]*temp3[j]*lamda[0][2];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp4[k]*lamda[0][3]+2*RTR11P1[k][j]*temp1[k]*temp4[j]*lamda[0][3];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp5[k]*lamda[0][4]+2*RTR11P1[k][j]*temp1[k]*temp5[j]*lamda[0][4];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp6[k]*lamda[0][5]+2*RTR11P1[k][j]*temp1[k]*temp6[j]*lamda[0][5];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp7[k]*lamda[0][6]+2*RTR11P1[k][j]*temp1[k]*temp7[j]*lamda[0][6];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp8[k]*lamda[0][7]+2*RTR11P1[k][j]*temp1[k]*temp8[j]*lamda[0][7];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp9[k]*lamda[0][8]+2*RTR11P1[k][j]*temp1[k]*temp9[j]*lamda[0][8];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp10[k]*lamda[0][9]+2*RTR11P1[k][j]*temp1[k]*temp10[j]*lamda[0][9];
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp11[k]*lamda[0][10]+2*RTR11P1[k][j]*temp1[k]*temp11[j]*lamda[0][10];
               
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp3[k]*lamda[1][2]+2*RTR11P1[k][j]*temp2[k]*temp3[j]*lamda[1][2];
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp4[k]*lamda[1][3]+2*RTR11P1[k][j]*temp2[k]*temp4[j]*lamda[1][3];
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp5[k]*lamda[1][4]+2*RTR11P1[k][j]*temp2[k]*temp5[j]*lamda[1][4];
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp6[k]*lamda[1][5]+2*RTR11P1[k][j]*temp2[k]*temp6[j]*lamda[1][5];
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp7[k]*lamda[1][6]+2*RTR11P1[k][j]*temp2[k]*temp7[j]*lamda[1][6];
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp8[k]*lamda[1][7]+2*RTR11P1[k][j]*temp2[k]*temp8[j]*lamda[1][7];
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp9[k]*lamda[1][8]+2*RTR11P1[k][j]*temp2[k]*temp9[j]*lamda[1][8];
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp10[k]*lamda[1][9]+2*RTR11P1[k][j]*temp2[k]*temp10[j]*lamda[1][9];
              func=func+ 2*RTR11P1[j][k]*temp2[j]*temp11[k]*lamda[1][10]+2*RTR11P1[k][j]*temp2[k]*temp11[j]*lamda[1][10];
           
              func=func+ 2*RTR11P1[j][k]*temp3[j]*temp4[k]*lamda[2][3]+ 2*RTR11P1[k][j]*temp3[k]*temp4[j]*lamda[2][3];
              func=func+ 2*RTR11P1[j][k]*temp3[j]*temp5[k]*lamda[2][4]+ 2*RTR11P1[k][j]*temp3[k]*temp5[j]*lamda[2][4];
              func=func+ 2*RTR11P1[j][k]*temp3[j]*temp6[k]*lamda[2][5]+ 2*RTR11P1[k][j]*temp3[k]*temp6[j]*lamda[2][5];
              func=func+ 2*RTR11P1[j][k]*temp3[j]*temp7[k]*lamda[2][6]+ 2*RTR11P1[k][j]*temp3[k]*temp7[j]*lamda[2][6];
              func=func+ 2*RTR11P1[j][k]*temp3[j]*temp8[k]*lamda[2][7]+ 2*RTR11P1[k][j]*temp3[k]*temp8[j]*lamda[2][7];
              func=func+ 2*RTR11P1[j][k]*temp3[j]*temp9[k]*lamda[2][8]+ 2*RTR11P1[k][j]*temp3[k]*temp9[j]*lamda[2][8];
              func=func+ 2*RTR11P1[j][k]*temp3[j]*temp11[k]*lamda[2][10]+ 2*RTR11P1[k][j]*temp3[k]*temp11[j]*lamda[2][10];
               
              func=func+ 2*RTR11P1[j][k]*temp4[j]*temp5[k]*lamda[3][4]+ 2*RTR11P1[k][j]*temp4[k]*temp5[j]*lamda[3][4];
              func=func+ 2*RTR11P1[j][k]*temp4[j]*temp6[k]*lamda[3][5]+ 2*RTR11P1[k][j]*temp4[k]*temp6[j]*lamda[3][5];
              func=func+ 2*RTR11P1[j][k]*temp4[j]*temp7[k]*lamda[3][6]+ 2*RTR11P1[k][j]*temp4[k]*temp7[j]*lamda[3][6];
              func=func+ 2*RTR11P1[j][k]*temp4[j]*temp8[k]*lamda[3][7]+ 2*RTR11P1[k][j]*temp4[k]*temp8[j]*lamda[3][7];
              func=func+ 2*RTR11P1[j][k]*temp4[j]*temp9[k]*lamda[3][8]+ 2*RTR11P1[k][j]*temp4[k]*temp9[j]*lamda[3][8];
              func=func+ 2*RTR11P1[j][k]*temp4[j]*temp10[k]*lamda[3][9]+ 2*RTR11P1[k][j]*temp4[k]*temp10[j]*lamda[3][9];
              func=func+ 2*RTR11P1[j][k]*temp4[j]*temp11[k]*lamda[3][10]+ 2*RTR11P1[k][j]*temp4[k]*temp11[j]*lamda[3][10];
            
              func=func+ 2*RTR11P1[j][k]*temp5[j]*temp6[k]*lamda[4][5]+ 2*RTR11P1[k][j]*temp5[k]*temp6[j]*lamda[4][5];
              func=func+ 2*RTR11P1[j][k]*temp5[j]*temp7[k]*lamda[4][6]+ 2*RTR11P1[k][j]*temp5[k]*temp7[j]*lamda[4][6];
              func=func+ 2*RTR11P1[j][k]*temp5[j]*temp8[k]*lamda[4][7]+ 2*RTR11P1[k][j]*temp5[k]*temp8[j]*lamda[4][7];
              func=func+ 2*RTR11P1[j][k]*temp5[j]*temp9[k]*lamda[4][8]+ 2*RTR11P1[k][j]*temp5[k]*temp9[j]*lamda[4][8];
              func=func+ 2*RTR11P1[j][k]*temp5[j]*temp10[k]*lamda[4][9]+ 2*RTR11P1[k][j]*temp5[k]*temp10[j]*lamda[4][9];
              func=func+ 2*RTR11P1[j][k]*temp5[j]*temp11[k]*lamda[4][10]+ 2*RTR11P1[k][j]*temp5[k]*temp11[j]*lamda[4][10];
             
            
              func=func+ 2*RTR11P1[j][k]*temp6[j]*temp7[k]*lamda[5][6]+ 2*RTR11P1[k][j]*temp6[k]*temp7[j]*lamda[5][6];
              func=func+ 2*RTR11P1[j][k]*temp6[j]*temp8[k]*lamda[5][7]+ 2*RTR11P1[k][j]*temp6[k]*temp8[j]*lamda[5][7];
              func=func+ 2*RTR11P1[j][k]*temp6[j]*temp9[k]*lamda[5][8]+ 2*RTR11P1[k][j]*temp6[k]*temp9[j]*lamda[5][8];
              func=func+ 2*RTR11P1[j][k]*temp6[j]*temp10[k]*lamda[5][9]+ 2*RTR11P1[k][j]*temp6[k]*temp10[j]*lamda[5][9];
              func=func+ 2*RTR11P1[j][k]*temp6[j]*temp11[k]*lamda[5][10]+ 2*RTR11P1[k][j]*temp6[k]*temp11[j]*lamda[5][10];
             
              func=func+ 2*RTR11P1[j][k]*temp7[j]*temp8[k]*lamda[6][7]+ 2*RTR11P1[k][j]*temp7[k]*temp8[j]*lamda[6][7];
              func=func+ 2*RTR11P1[j][k]*temp7[j]*temp9[k]*lamda[6][8]+ 2*RTR11P1[k][j]*temp7[k]*temp9[j]*lamda[6][8];
              func=func+ 2*RTR11P1[j][k]*temp7[j]*temp10[k]*lamda[6][9]+ 2*RTR11P1[k][j]*temp7[k]*temp10[j]*lamda[6][9];
              func=func+ 2*RTR11P1[j][k]*temp7[j]*temp11[k]*lamda[6][10]+ 2*RTR11P1[k][j]*temp7[k]*temp11[j]*lamda[6][10];
           
              func=func+ 2*RTR11P1[j][k]*temp8[j]*temp9[k]*lamda[7][8]+ 2*RTR11P1[k][j]*temp8[k]*temp9[j]*lamda[7][8];
              func=func+ 2*RTR11P1[j][k]*temp8[j]*temp10[k]*lamda[7][9]+ 2*RTR11P1[k][j]*temp8[k]*temp10[j]*lamda[7][9];
              func=func+ 2*RTR11P1[j][k]*temp8[j]*temp11[k]*lamda[7][10]+ 2*RTR11P1[k][j]*temp8[k]*temp11[j]*lamda[7][10];
          
              func=func+ 2*RTR11P1[j][k]*temp9[j]*temp10[k]*lamda[8][9]+ 2*RTR11P1[k][j]*temp8[k]*temp10[j]*lamda[8][9];
              func=func+ 2*RTR11P1[j][k]*temp9[j]*temp11[k]*lamda[8][10]+ 2*RTR11P1[k][j]*temp8[k]*temp11[j]*lamda[8][10];
               
              func=func+ 2*RTR11P1[j][k]*temp10[j]*temp11[k]*lamda[9][10]+ 2*RTR11P1[k][j]*temp10[k]*temp11[j]*lamda[9][10];
  
            }
    }
    sum=sum-0.5*func;
   
    Free(temp1);
    Free(temp2);
    Free(temp3);
    Free(temp4);
    Free(temp5);
    Free(temp6);
    Free(temp7);
    Free(temp8);
    Free(temp9);
    Free(temp10);
    Free(temp11);
    return (sum);
}

double target_sig (double b0, double b1, double b2, double b3, double b4, double b5, double b6, double b7,
                   double b8, double b9, double b10, double b11,double b12, double b14, double b15,double b16,double b17, double b18,
                   double b19,double b20, double b21, double b22, double b23,double *phi1, double *phi2, double *phi3, double *phi4, double *phi5, 
                   double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11,  double x, int *Cores, 
                   double *year1, double *year2,double *year3,double *year4, double *year5, double *year6, double *year7, double *year8, double *year9,
                   double *year10, double *x1, double *x2, double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9,
                    double *x10, double *x11, double *x12, double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
     a=1/(x);
  
    for (int i = 0; i < N; i++)
    {
         mu=b0 + b1*x1[i]+b2*x2[i]+b14*x3[i]+b15*x4[i]+b16*x5[i]+b17*x6[i]+b18*x7[i]+b19*x8[i]+b20*x9[i]+b21*x10[i]+b22*x11[i]+b23*x12[i]
           +b3*year1[i]+b4*year2[i]+b5*year3[i]+b6*year4[i]+b7*year5[i]+b8*year6[i]+b9*year7[i]+b10*year8[i]+b11*year9[i]+b12*year10[i]
           +phi1[Cores[i]]+phi2[Cores[i]]*year1[i]+phi3[Cores[i]]*year2[i]+phi4[Cores[i]]*year3[i]+phi5[Cores[i]]*year4[i]+phi6[Cores[i]]*year5[i]
           +phi7[Cores[i]]*year6[i]+phi8[Cores[i]]*year7[i]+phi9[Cores[i]]*year8[i]+phi10[Cores[i]]*year9[i]+phi11[Cores[i]]*year10[i];
     
        lam=exp(-mu/x);
        sum=sum+s[i]*log(lam*a*pow(y[i],(a-1)))-lam*pow(y[i],a);
    }
    
    sum=sum+dgamma(x,0.001,1000,1);
    
    return sum;
}


/***** for lamda matrix ***********************/
void target_lamda(double **matrix, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, 
                   double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11,  int no_regions, double **RTR11P1, int p)
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
            temp[0][2]=temp[0][2]+ phi1[i]*phi3[j]*RTR11P1[i][j];
            temp[0][3]=temp[0][3]+ phi1[i]*phi4[j]*RTR11P1[i][j];
            temp[0][4]=temp[0][4]+ phi1[i]*phi5[j]*RTR11P1[i][j];
            temp[0][5]=temp[0][5]+ phi1[i]*phi6[j]*RTR11P1[i][j];
            temp[0][6]=temp[0][6]+ phi1[i]*phi7[j]*RTR11P1[i][j];
            temp[0][7]=temp[0][7]+ phi1[i]*phi8[j]*RTR11P1[i][j];
            temp[0][8]=temp[0][8]+ phi1[i]*phi9[j]*RTR11P1[i][j];
            temp[0][9]=temp[0][9]+ phi1[i]*phi10[j]*RTR11P1[i][j];
            temp[0][10]=temp[0][10]+ phi1[i]*phi11[j]*RTR11P1[i][j];
            
            temp[1][1]=temp[1][1]+ phi2[i]*phi2[j]*RTR11P1[i][j];
            temp[1][2]=temp[1][2]+ phi2[i]*phi3[j]*RTR11P1[i][j];
            temp[1][3]=temp[1][3]+ phi2[i]*phi4[j]*RTR11P1[i][j];
            temp[1][4]=temp[1][4]+ phi2[i]*phi5[j]*RTR11P1[i][j];
            temp[1][5]=temp[1][5]+ phi2[i]*phi6[j]*RTR11P1[i][j];
            temp[1][6]=temp[1][6]+ phi2[i]*phi7[j]*RTR11P1[i][j];
            temp[1][7]=temp[1][7]+ phi2[i]*phi8[j]*RTR11P1[i][j];
            temp[1][8]=temp[1][8]+ phi2[i]*phi9[j]*RTR11P1[i][j];
            temp[1][9]=temp[1][9]+ phi2[i]*phi10[j]*RTR11P1[i][j];
            temp[1][10]=temp[1][10]+ phi2[i]*phi11[j]*RTR11P1[i][j];
          

            temp[2][2]=temp[2][2]+ phi3[i]*phi3[j]*RTR11P1[i][j];
            temp[2][3]=temp[2][3]+ phi3[i]*phi4[j]*RTR11P1[i][j];
            temp[2][4]=temp[2][4]+ phi3[i]*phi5[j]*RTR11P1[i][j];
            temp[2][5]=temp[2][5]+ phi3[i]*phi6[j]*RTR11P1[i][j];
            temp[2][6]=temp[2][6]+ phi3[i]*phi7[j]*RTR11P1[i][j];
            temp[2][7]=temp[2][7]+ phi3[i]*phi8[j]*RTR11P1[i][j];
            temp[2][8]=temp[2][8]+ phi3[i]*phi9[j]*RTR11P1[i][j];
            temp[2][9]=temp[2][9]+ phi3[i]*phi10[j]*RTR11P1[i][j];
            temp[2][10]=temp[2][10]+ phi3[i]*phi11[j]*RTR11P1[i][j];
          

            temp[3][3]=temp[3][3]+ phi4[i]*phi4[j]*RTR11P1[i][j];
            temp[3][4]=temp[3][4]+ phi4[i]*phi5[j]*RTR11P1[i][j];
            temp[3][5]=temp[3][5]+ phi4[i]*phi6[j]*RTR11P1[i][j];
            temp[3][6]=temp[3][6]+ phi4[i]*phi7[j]*RTR11P1[i][j];
            temp[3][7]=temp[3][7]+ phi4[i]*phi8[j]*RTR11P1[i][j];
            temp[3][8]=temp[3][8]+ phi4[i]*phi9[j]*RTR11P1[i][j];
            temp[3][9]=temp[3][9]+ phi4[i]*phi10[j]*RTR11P1[i][j];
            temp[3][10]=temp[3][10]+ phi4[i]*phi11[j]*RTR11P1[i][j];
          
            temp[4][4]=temp[4][4]+ phi5[i]*phi5[j]*RTR11P1[i][j];
            temp[4][5]=temp[4][5]+ phi5[i]*phi6[j]*RTR11P1[i][j];
            temp[4][6]=temp[4][6]+ phi5[i]*phi7[j]*RTR11P1[i][j];
            temp[4][7]=temp[4][7]+ phi5[i]*phi8[j]*RTR11P1[i][j];
            temp[4][8]=temp[4][8]+ phi5[i]*phi9[j]*RTR11P1[i][j];
            temp[4][9]=temp[4][9]+ phi5[i]*phi10[j]*RTR11P1[i][j];
            temp[4][10]=temp[4][10]+ phi5[i]*phi11[j]*RTR11P1[i][j];
  
            temp[5][5]=temp[5][5]+ phi6[i]*phi6[j]*RTR11P1[i][j];
            temp[5][6]=temp[5][6]+ phi6[i]*phi7[j]*RTR11P1[i][j];
            temp[5][7]=temp[5][7]+ phi6[i]*phi8[j]*RTR11P1[i][j];
            temp[5][8]=temp[5][8]+ phi6[i]*phi9[j]*RTR11P1[i][j];
            temp[5][9]=temp[5][9]+ phi6[i]*phi10[j]*RTR11P1[i][j];
            temp[5][10]=temp[5][10]+ phi6[i]*phi11[j]*RTR11P1[i][j];
 
            temp[6][6]=temp[6][6]+ phi7[i]*phi7[j]*RTR11P1[i][j];
            temp[6][7]=temp[6][7]+ phi7[i]*phi8[j]*RTR11P1[i][j];
            temp[6][8]=temp[6][8]+ phi7[i]*phi9[j]*RTR11P1[i][j];
            temp[6][9]=temp[6][9]+ phi7[i]*phi10[j]*RTR11P1[i][j];
            temp[6][10]=temp[6][10]+ phi7[i]*phi11[j]*RTR11P1[i][j];
 
            temp[7][7]=temp[7][7]+ phi8[i]*phi8[j]*RTR11P1[i][j];
            temp[7][8]=temp[7][8]+ phi8[i]*phi9[j]*RTR11P1[i][j];
            temp[7][9]=temp[7][9]+ phi8[i]*phi10[j]*RTR11P1[i][j];
            temp[7][10]=temp[7][10]+ phi8[i]*phi11[j]*RTR11P1[i][j];
 
            temp[8][8]=temp[8][8]+ phi9[i]*phi9[j]*RTR11P1[i][j];
            temp[8][9]=temp[8][9]+ phi9[i]*phi10[j]*RTR11P1[i][j];
            temp[8][10]=temp[8][10]+ phi9[i]*phi11[j]*RTR11P1[i][j];
  
            temp[9][9]=temp[9][9]+ phi10[i]*phi10[j]*RTR11P1[i][j];
            temp[9][10]=temp[9][10]+ phi10[i]*phi11[j]*RTR11P1[i][j];
 
            temp[10][10]=temp[10][10]+ phi11[i]*phi11[j]*RTR11P1[i][j];
             
           
        }
        
    }
    
   for (i=0;i<p;i++){
       for (j=0;j<i;j++){
           temp[i][j]=temp[j][i];
       }
   }


    transpose(matrix, temp, p, p);   
    freematd(temp,p, p);
    
}


/********* Metroplis-Hasting for rho and beta0  ***********************/

double metroplis_b0(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,
                   double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,
                    double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, 
                   double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11,  double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7,
                    double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);

    y1=b0+0.1*rnorm(0,1);
  
    temp1=target_b(y1,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,  sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b0;}
    }
    else {newx=b0;}
    return (newx);
}


double metroplis_b1(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11,  double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,   double *x9, double *x10, double *x11, double *x12,double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b1+0.1*rnorm(0,1);
    
    temp1=target_b(b0,y1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b1;}
    }
    else {newx=b1;}
    return (newx);
}


double metroplis_b2(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2, double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b2+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,y1,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b2;}
    }
    else {newx=b2;}
    return (newx);
}

double metroplis_b3(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12,  double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b3+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,y1,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b3;}
    }
    else {newx=b3;}
    return (newx);
}

double metroplis_b4(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b4+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,y1,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b4;}
    }
    else {newx=b4;}
    return (newx);
}

double metroplis_b5(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b5+0.1*rnorm(0,1);

    temp1=target_b(b0,b1,b2,b3,b4,y1,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b5;}
    }
    else {newx=b5;}
    return (newx);
}

double metroplis_b6(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b6+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,y1,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b6;}
    }
    else {newx=b6;}
    return (newx);
}

double metroplis_b7(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b7+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,y1,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b7;}
    }
    else {newx=b7;}
    return (newx);
}

double metroplis_b8(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b8+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,y1,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b8;}
    }
    else {newx=b8;}
    return (newx);
}

double metroplis_b9(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,  double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b9+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,y1,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b9;}
    }
    else {newx=b9;}
    return (newx);
}

double metroplis_b10(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b10+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,y1,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b9;}
    }
    else {newx=b9;}
    return (newx);
}


double metroplis_b11(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,  double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b11+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,y1,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b11;}
    }
    else {newx=b11;}
    return (newx);
}

double metroplis_b12(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,  double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);   
    y1=b12+0.1*rnorm(0,1);
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,y1,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b12;}
    }
    else {newx=b12;}
    return (newx);
}


double metroplis_b14(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,  double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b14+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,y1,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b14;}
    }
    else {newx=b14;}
    return (newx);
}

double metroplis_b15(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,  double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b15+0.01*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,y1,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b15;}
    }
    else {newx=b15;}
    return (newx);
}

double metroplis_b16(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12,  double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b16+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,y1,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b16;}
    }
    else {newx=b16;}
    return (newx);
}

double metroplis_b17(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b17+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,y1,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b17;}
    }
    else {newx=b17;}
    return (newx);
}

double metroplis_b18(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b18+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,y1,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b18;}
    }
    else {newx=b18;}
    return (newx);
}

double metroplis_b19(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b19+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,y1,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b19;}
    }
    else {newx=b19;}
    return (newx);
}


double metroplis_b20(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b20+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,y1,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b20;}
    }
    else {newx=b20;}
    return (newx);
}

double metroplis_b21(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b21+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,y1,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b21;}
    }
    else {newx=b21;}
    return (newx);
}

double metroplis_b22(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b22+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,y1,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b22;}
    }
    else {newx=b22;}
    return (newx);
}

double metroplis_b23(double b0, double b1,double b2,double b3, double b4, double b5, double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11, double sig,  int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y, int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);
    
    y1=b23+0.1*rnorm(0,1);
    
    temp1=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,y1,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_b(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi1, phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=b23;}
    }
    else {newx=b23;}
    return (newx);
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

void metroplis_phi(double b0, double b1,double b2,double b3, double b4,double b5,double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23,double *phi_curr1, double *phi_curr2, double *phi_curr3, double *phi_curr4, double *phi_curr5, double *phi_curr6, double *phi_curr7, double *phi_curr8, double *phi_curr9, double *phi_curr10,double *phi_curr11, double sig, double **lamda, int *Cores, double *year1, double *year2,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *year3, double *year4, double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,   double *y,int *s, int N, int no_regions, double **RTR11P1,double x1_prev, double x2_prev,double x3_prev, double x4_prev,double x5_prev, double x6_prev, double x7_prev, double x8_prev, double x9_prev, double x10_prev, double x11_prev, int j)
{
    double lgratio=0, u ,temp1=0, temp2=0;
    double x1_curr, x2_curr, x3_curr, x4_curr, x5_curr, x6_curr, x7_curr, x8_curr, x9_curr, x10_curr, x11_curr;
   
    
    u=runif(0,1);
   
    x1_curr=x1_prev+rnorm(0,0.1);
    x2_curr=x2_prev+rnorm(0,0.1);
    x3_curr=x3_prev+rnorm(0,0.1);
    x4_curr=x4_prev+rnorm(0,0.1);
    x5_curr=x5_prev+rnorm(0,0.1);
    x6_curr=x6_prev+rnorm(0,0.1);
    x7_curr=x7_prev+rnorm(0,0.1);
    x8_curr=x8_prev+rnorm(0,0.1);
    x9_curr=x9_prev+rnorm(0,0.1);
    x10_curr=x10_prev+rnorm(0,0.1);
    x11_curr=x11_prev+rnorm(0,0.1);

    
    temp1=target_phi(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi_curr1,phi_curr2,phi_curr3,phi_curr4,phi_curr5,phi_curr6,phi_curr7,phi_curr8,phi_curr9,
                      phi_curr10,phi_curr11,sig, Cores, year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, 
                      x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N, no_regions,RTR11P1,lamda, x1_curr, x2_curr,x3_curr,x4_curr,x5_curr,x6_curr,x7_curr,x8_curr,x9_curr,x10_curr,x11_curr,j);
    temp2=target_phi(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi_curr1,phi_curr2,phi_curr3,phi_curr4,phi_curr5,phi_curr6,phi_curr7,phi_curr8,phi_curr9,
                      phi_curr10,phi_curr11,
                      sig, Cores, year1,year2,year3,year4, year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N, no_regions,RTR11P1,lamda, 
                       x1_prev,x2_prev,x3_prev,x4_prev,x5_prev,x6_prev,x7_prev,x8_prev,x9_prev,x10_prev,x11_prev,j);
   lgratio=temp1-temp2;
    // if(j==0) printf("lgratio is %f \n", lgratio);
    if (u>0)
    { if (lgratio>=log(u))
    {   
        phi_curr1[j]=x1_curr;
        phi_curr2[j]=x2_curr;
        phi_curr3[j]=x3_curr;
        phi_curr4[j]=x4_curr;
        phi_curr5[j]=x5_curr;
        phi_curr6[j]=x6_curr;
        phi_curr7[j]=x7_curr;
        phi_curr8[j]=x8_curr;
        phi_curr9[j]=x9_curr;
        phi_curr10[j]=x10_curr;
        phi_curr11[j]=x11_curr;
    }
    else
    {
        phi_curr1[j]=x1_prev;
        phi_curr2[j]=x2_prev;
        phi_curr3[j]=x3_prev;
        phi_curr4[j]=x4_prev;
        phi_curr5[j]=x5_prev;
        phi_curr6[j]=x6_prev;
        phi_curr7[j]=x7_prev;
        phi_curr8[j]=x8_prev;
        phi_curr9[j]=x9_prev;
        phi_curr10[j]=x10_prev;
        phi_curr11[j]=x11_prev;

        
    }
    }
    else
        
    {   phi_curr1[j]=x1_prev;
        phi_curr2[j]=x2_prev;
        phi_curr3[j]=x3_prev;
        phi_curr4[j]=x4_prev;
        phi_curr5[j]=x5_prev;
        phi_curr6[j]=x6_prev;
        phi_curr7[j]=x7_prev;
        phi_curr8[j]=x8_prev;
        phi_curr9[j]=x9_prev;
        phi_curr10[j]=x10_prev;
        phi_curr11[j]=x11_prev;

    }
    
}



double metroplis_sig( double b0, double b1, double b2, double b3,double b4,double b5,double b6,double b7, double b8, double b9, double b10, double b11,double b12,double b14, double b15,double b16,double b17,double b18,double b19, double b20, double b21, double b22, double b23,double *phi_curr1, double *phi_curr2,double *phi_curr3, double *phi_curr4, double *phi_curr5, double *phi_curr6, double *phi_curr7, double *phi_curr8, double *phi_curr9, double *phi_curr10,double *phi_curr11,double x_prev, int *Cores, double *year1, double *year2,double *year3, double *year4,double *year5, double *year6, double *year7, double *year8, double *year9,double *year10, double *x1, double *x2, double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *x9, double *x10, double *x11, double *x12, double *y,int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+rnorm(0,0.1);
    if (y_new>0){
    temp1=target_sig(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi_curr1, phi_curr2,phi_curr3,phi_curr4,phi_curr5,phi_curr6,phi_curr7,phi_curr8,phi_curr9,phi_curr10,phi_curr11, y_new,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    temp2=target_sig(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,phi_curr1,phi_curr2,phi_curr3,phi_curr4,phi_curr5,phi_curr6,phi_curr7,phi_curr8,phi_curr9,phi_curr10,phi_curr11, x_prev,Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    
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


/*double LL (double b0, double b1, double b2, double b3,double b4,double b5,double b6,double b7, double b8, double b9, double b10, double b11,double b12,double  double b14, double b15,double b16,double b17,double b18,double b19,double b20, double b21, double b22, double b23,double *phi1, double *phi2,double sig, int *Cores, double *year,double *x1, double *x2,double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8,  double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
    for (int i = 0; i < N; i++){
        mu=b0 + b1*x1[i]+b2*x2[i]+b3*year1[i]+b4*year2[i]+b5*year3[i]+b6*year4[i]+b7*year5[i]+b8*year6[i]+b9*year7[i]+b10*year8[i]+b11*year9[i]+b12*year10[i]+b13*year11[i]+b14*x3[i]+ b15*x4[i]+b16*x5[i]+b17*x6[i]+b18*x7[i]+b19*x8[i]+phi1[Cores[i]]+phi2[Cores[i]]*(year1[i]+year2[i]+year3[i]+year4[i]+year5[i]+year6[i]+year7[i]+year8[i]+year9[i]+year10[i]+year11[i]);
        lam=exp(-mu/sig);
        a=1/(sig);
        sum=sum+s[i]*log(lam*a*pow(y[i],(a-1)))-lam*pow(y[i],a);
        
    }
    return -2*sum;
}*/


void weibul(double *b0, double *b1, double *b2, double *b3,double *b4,double *b5,double *b6,double *b7, double *b8, double *b9, double *b10, double *b11,double *b12, double *b14, double *b15,double *b16,double *b17,double *b18,double *b19,double *b20, double *b21, double *b22, double *b23,double *sig, double *lamda_p0, double *phi1, double *phi2,double *phi3, double *phi4, double *phi5, double *phi6, double *phi7, double *phi8, double *phi9, double *phi10,double *phi11,  double *vecD, double *vecC,double *y, int *s, int *N_p, double *x1, double *x2, double  *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, double *x11, double *x12,  double *year,int *Cores, int *nend, int *no_regions_p,int *dims,double *lamda_prec,double *sumLnF,double *phi_sample,int *n_burn){
	int i, j;
    int p=dims[0];
    int N=N_p[0];
    int no_regions=no_regions_p[0];
    double **temp_lamda, **tempsum, *tempsum_INV, **lamda, **lamda_prec_mat,**D,**C,**RTR11P1;
    double *year1, *year2, *year3, *year4, *year5, *year6, *year7, *year8,*year9,*year10;
    double x1_prev,x2_prev,x3_prev,x4_prev,x5_prev,x6_prev, x7_prev, x8_prev, x9_prev,x10_prev, x11_prev;
    
    year1=create_vectord(N);
    year2=create_vectord(N);
    year3=create_vectord(N);
    year4=create_vectord(N);
    year5=create_vectord(N);
    year6=create_vectord(N);
    year7=create_vectord(N);
    year8=create_vectord(N);
    year9=create_vectord(N);
    year10=create_vectord(N);

  
  for(int i =0;i<N;i++){
        year1[i]=(year[i]==6);
        year2[i]=(year[i]==7);
        year3[i]=(year[i]==8);
        year4[i]=(year[i]==9);
       year5[i]=(year[i]==10);
       year6[i]=(year[i]==11);
       year7[i]=(year[i]==12);
       year8[i]=(year[i]==13);
       year9[i]=(year[i]==14);
       year10[i]=(year[i]==15);
    }
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
    
	{

     b0[i]=metroplis_b0(b0[i-1],b1[i-1],b2[i-1],b3[i-1],b4[i-1],b5[i-1],b6[i-1],b7[i-1],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1, year2,year3,year4,year5,year6,year7,year8,year9,year10,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b1[i]=metroplis_b1(b0[i],b1[i-1],b2[i-1],b3[i-1],b4[i-1],b5[i-1],b6[i-1],b7[i-1],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2, phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
     b2[i]=metroplis_b2(b0[i],b1[i],b2[i-1],b3[i-1],b4[i-1],b5[i-1],b6[i-1],b7[i-1],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b3[i]=metroplis_b3(b0[i],b1[i],b2[i],b3[i-1],b4[i-1],b5[i-1],b6[i-1],b7[i-1],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    b4[i]=metroplis_b4(b0[i],b1[i],b2[i],b3[i],b4[i-1],b5[i-1],b6[i-1],b7[i-1],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
     b5[i]=metroplis_b5(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i-1],b6[i-1],b7[i-1],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
     b6[i]=metroplis_b6(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i-1],b7[i-1],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
     b7[i]=metroplis_b7(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i-1],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
    b8[i]=metroplis_b8(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i-1],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
    b9[i]=metroplis_b9(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i-1],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
    b10[i]=metroplis_b10(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i-1],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
    b11[i]=metroplis_b11(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i-1],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);

   b12[i]=metroplis_b12(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i-1],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
        
    b14[i]=metroplis_b14(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i-1],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
    b15[i]=metroplis_b15(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i-1],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
    b16[i]=metroplis_b16(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i-1],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);

    b17[i]=metroplis_b17(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i-1],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
        
    b18[i]=metroplis_b18(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i-1],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
      
     b19[i]=metroplis_b19(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i-1],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);

    b20[i]=metroplis_b20(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i-1],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);

    b21[i]=metroplis_b21(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i],b21[i-1],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
   
    b22[i]=metroplis_b22(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i],b21[i],b22[i-1],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    
    b23[i]=metroplis_b23(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i],b21[i],b22[i],b23[i-1],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);
    
    for (j = 0; j < no_regions; j++)
	    {
            x1_prev=phi1[j];
            x2_prev=phi2[j];
            x3_prev=phi3[j];
            x4_prev=phi4[j];
            x5_prev=phi5[j];
            x6_prev=phi6[j];
            x7_prev=phi7[j];
            x8_prev=phi8[j];
            x9_prev=phi9[j];
            x10_prev=phi10[j];
            x11_prev=phi11[j];
            
            metroplis_phi(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i],b21[i],b22[i],b23[i],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],lamda, Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N,no_regions,RTR11P1,x1_prev,x2_prev,x3_prev,x4_prev, x5_prev, x6_prev,x7_prev,x8_prev,x9_prev,x10_prev, x11_prev, j);
        }
    //printf("Before center %d is %f %f %f \n",i+1,phi1[0],phi1[1],phi1[2]);
    //printf("Before center %d is %f %f %f \n",i+1,phi2[0],phi2[1],phi2[2]);
      phi1=center(phi1,no_regions);
      phi2=center(phi2,no_regions);
      phi3=center(phi3,no_regions);
      phi4=center(phi4,no_regions);
      phi5=center(phi5,no_regions);
      phi6=center(phi6,no_regions);
      phi7=center(phi7,no_regions);
      phi8=center(phi8,no_regions);
      phi9=center(phi9,no_regions);
      phi10=center(phi10,no_regions);
      phi11=center(phi11,no_regions);
     
      
        
      sig[i]=metroplis_sig(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i],b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],b20[i],b21[i],b22[i],b23[i],phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, sig[i-1],Cores,year1,year2,year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y,s,N);

      target_lamda(temp_lamda,phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11, no_regions,RTR11P1,p);
      matsum (tempsum, lamda_prec_mat, temp_lamda, p, p);
      tempsum_INV=inverse(tempsum_INV,p,tempsum);
      rwishart(lamda_p0,p+no_regions,tempsum_INV,dims);
      vectormat(lamda,lamda_p0,p,p);
      //if(i>=n_burn[0]) sumLnF[0]=sumLnF[0]+LL(b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i], b7[i],b8[i],b9[i],b10[i],b11[i],b12[i],b14[i],b15[i],b16[i],b17[i],b18[i],b19[i],phi1,phi2, sig[i],Cores,year1,year2, year3,year4,year5,year6,year7,year8,year9,year10, x1,x2,x3,x4,x5,x6,x7,x8,y,s,N);
     for(int k=0;k<no_regions;k++){
          phi_sample[k*nend[0]+i]=phi1[k];
          phi_sample[(k+no_regions)*nend[0]+i]=phi2[k];
          phi_sample[(k+2*no_regions)*nend[0]+i]=phi3[k];
          phi_sample[(k+3*no_regions)*nend[0]+i]=phi4[k];
          phi_sample[(k+4*no_regions)*nend[0]+i]=phi5[k];
          phi_sample[(k+5*no_regions)*nend[0]+i]=phi6[k];
          phi_sample[(k+6*no_regions)*nend[0]+i]=phi7[k];
          phi_sample[(k+7*no_regions)*nend[0]+i]=phi8[k];
          phi_sample[(k+8*no_regions)*nend[0]+i]=phi9[k];
          phi_sample[(k+9*no_regions)*nend[0]+i]=phi10[k];
          phi_sample[(k+10*no_regions)*nend[0]+i]=phi11[k];
      

        }
    }
    
 Free(year1);
 Free(year2);
 Free(year3);
 Free(year4);
 Free(year5);
 Free(year6);
 Free(year7);
 Free(year8);
 Free(year9);
 Free(year10);
 
    freematd(tempsum,p, p);
    Free(tempsum_INV);
   
    freematd(lamda_prec_mat, p, p);
    freematd(RTR11P1,no_regions,no_regions);
    freematd(C, no_regions, no_regions);
    freematd(D, no_regions, no_regions);
    freematd(lamda, p,p);
    PutRNGstate( );

}






