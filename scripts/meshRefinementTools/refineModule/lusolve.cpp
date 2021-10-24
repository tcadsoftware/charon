
#include <cmath>
#include <iostream>

//Very simple lu decomposition with partial pivoting


int ludcmp(double a[][4], int n, int indx[], int *d) 

{ 

  double *vv,big,sum,dum,small;
  int imax=0;

  small=1.0e-16; 

  vv = new double[100];

  *d = 1;

  for(int i=0 ; i<n ; ++i) 
    { 
      big = 0; 
      for(int j=0 ; j<n ; j++) 
	if(fabs(a[i][j]) > big) 
	  big = fabs(a[i][j]);
      if(big == 0)
	{
	  std::cout<<" MATRIX HAS ZERO ROW   "<<i<<std::endl;
	  for(int j=0 ; j<n ; ++j)
	    std::cout<<fabs(a[i][j])<<"    ";
	  std::cout<<std::endl;
	  return(0);  //Zero row 
	}
      vv[i] = 1/big; 
    } 
  
  for(int j=0 ; j<n ; j++) 
    { 
      if(j != 0) 
	for(int i=0 ; i<j ; i++) 
	  { 
	    sum=a[i][j]; 
	    if(i != 0) 
	      for(int k=0 ; k<i ; k++) 
		sum-=(a[i][k]*a[k][j]); 
	    a[i][j]=sum; 
	  } 
      big = 0; 
      for(int i=j ; i<n ; i++) 
	{ 
	  sum=a[i][j]; 
	  if(j != 0) 
	    for(int k=0 ; k<j ; k++) 
	      sum-=a[i][k]*a[k][j]; 
	  a[i][j]=sum; 
	  dum=vv[i]*fabs(sum); 
	  if(dum >= big) 
	    { 
	      big=dum; imax=i; 
	    } 
	} 
      if(j != imax)
	{ 
	  for(int k=0 ; k<n ; k++) 
	    { 
	      dum=a[imax][k]; 
	      a[imax][k]=a[j][k]; 
	      a[j][k]=dum; 
	    } 
	  *d=-*d;
	  dum=vv[imax]; 
	  vv[imax]=vv[j];
	  vv[j]=dum; 
	} 
      indx[j]=imax;
      if(a[j][j] == 0) 
	a[j][j]=small;
      
      if(j != n-1) 
	{ 
	  dum=1/a[j][j]; 
	  for(int i=j+1 ; i<n ; ++i) 
	    a[i][j]=a[i][j]*dum; 
	} 
    }
  
  delete [] vv;

  return 1; 
} 


void lubksb( double a[][4], int n, int indx[], double b[]) 
{ 

  double sum;

  int ii = -1; 

  for(int i=0 ; i<n ; ++i) 
    { 
      int ip = indx[i]; 
      sum = b[ip]; 
      b[ip] = b[i]; 
      if(ii != -1) 
	for(int j=ii ; j<i ; ++j) 
	  sum=sum-a[i][j]*b[j]; 
      else
	if(sum!=0.0) 
	  ii=i; 
      b[i]=sum; 
    } 

  for(int i=n-1 ; i>=0 ; --i) 
    { 
      sum = b[i]; 
      if(i != n-1) 
	for(int j=i+1 ; j<n ; ++j) 
	  sum=sum-a[i][j]*b[j]; 
      b[i] = sum/a[i][i]; 
    } 
  return; 
} 

 //************ ** lusolve ** ************* 
 // Solve a linear set of equations: A x = b 
 // Original matrix A will be destroyed by this operation. 
 // Returns 0 if matrix is singular, 1 otherwise.

int lusolve(double a[][4], int n, double b[]) 
{ 
  
  int *indx; 
  int d; 
  
  indx = new int[n];

  if(ludcmp(a,n,indx,&d) == 0) 
    return 0; 
  
  lubksb(a,n,indx,b); 

  delete [] indx;

  return 1; 

} 
