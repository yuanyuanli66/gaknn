#include <malloc.h>
#include <math.h>
#include "ga_knn.h"


#define SQR(a) ((a)*(a));
void crank(int n, float w[], float *s);
void sort2(int n, float arr[], float brr[]);

float pearson_correlation(float *x,float *y,int n) {

   register int l;
   float sumXY=0;
   float sumX=0; 
   float sumY=0;
   float sumX2=0;
   float sumY2=0;
   for (l=0; l<n; l++) {
      sumXY +=x[l]*y[l];
      sumX  +=x[l];
      sumY  +=y[l];
      sumX2 +=x[l]*x[l];
      sumY2 +=y[l]*y[l];
   }
   float r=(n*sumXY-sumX*sumY)/sqrt((n*sumX2-sumX*sumX)*(n*sumY2-sumY*sumY));
   if (fabs(r)>1) r=0;
   return (r); 
}

float spearman_correlation(float *x,float *y,int n) {

   register int j;
   float sg,sf,fac,en3n,en;
   float x2[n+1],y2[n+1];

   for (j=0; j<n; j++)  {
      y2[j+1]=y[j];
      x2[j+1]=x[j];
   } 
   sort2(n,x2,y2);
   crank(n,x2,&sf);
   sort2(n,y2,x2);
   crank(n,y2,&sg);

   float d=0.0; for (j=1; j<=n; j++) d += SQR(x2[j]-y2[j]);

   en=n;
   en3n=en*en*en-en;
   fac=(1.0-sf/en3n)*(1.0-sg/en3n);
   float r=(1.0-(6.0/en3n)*(d+(sf+sg)/12.0))/sqrt(fac);
   return (r);
}

