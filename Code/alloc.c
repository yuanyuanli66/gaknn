
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include "ga_knn.h"

Similarity *alloc_sim(int num) {

   Similarity *tmp=(Similarity *)malloc(num*sizeof(Similarity));
   if (!tmp) { printf("distance matrix malloc failed!\n"); exit(1); }
   return (tmp);
}

Fitness **alloc_fitness_fitness(int size1,int size2) {

   register int i;
   Fitness **tmp=(Fitness **)calloc(size1,sizeof(Fitness *));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }

   tmp[0]=(Fitness *)calloc(size1*size2,sizeof(Fitness));
   if (tmp[0]==0) { printf("bit calloc failed!\n"); exit(1); }

   /* set up vector pointers */
   for (i=1; i<size1; i++) {
      tmp[i]=tmp[0] + (size2 * i);
   }
   return (tmp);
}

Fitness *alloc_fitness (int pop_size) {

   Fitness *tmp=(Fitness *)malloc(pop_size*sizeof(Fitness));
   if (!tmp) {
      printf("malloc failed for fitness !\n"); exit(1);
   }
   return (tmp);
}

Wheel *alloc_wheel (int pop_size) {

   Wheel *tmp=(Wheel *)malloc(pop_size*sizeof(Wheel));
   if (!tmp) {
      printf("malloc failed for fitness !\n"); exit(1);
   }
   return (tmp);
}

CumulativePred *alloc_cum (int size1) {

   CumulativePred *tmp=(CumulativePred *)malloc(size1*sizeof(CumulativePred));
   if (!tmp) {
      printf("malloc failed for fitness !\n"); exit(1);
   }
   return (tmp);
}

float *alloc_float (int pop_size) {

   float *tmp=(float *)malloc(pop_size*sizeof(float));
   if (!tmp) { printf("malloc failed for float!\n"); exit(1); }
   return (tmp);
}

int *alloc_int(int size1) {

   int *tmp=(int *)calloc(size1,sizeof(int));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }
   return (tmp);
}

int **alloc_int_int(int size1,int size2) {

   register int i;
   int **tmp=(int **)malloc(size1*sizeof(int *));
   if (!tmp)   { printf("tmp malloc failed!\n"); exit(1); }

   tmp[0]=(int *)calloc(size1*size2,sizeof(int));
   if (tmp[0]==0) { printf("bit calloc failed!\n"); exit(1); }

   /* set up vector pointers */
   for (i=1; i<size1; i++) {
      tmp[i]=tmp[0] + (size2 * i);
   }
   return (tmp);
}

float **alloc_float_float(int size1,int size2) {

   register int i;
   float **tmp=(float **)malloc(size1*sizeof(float *));
   if (!tmp)   { printf("tmp malloc failed!\n"); exit(1); }

   tmp[0]=(float *)calloc(size1*size2,sizeof(float));
   if (tmp[0]==0) { printf("bit calloc failed!\n"); exit(1); }

   /* set up vector pointers */
   for (i=1; i<size1; i++) {
      tmp[i]=tmp[0] + (size2 * i);
   }
   return (tmp);
}

char *alloc_char(int size1) {

   char *tmp=(char *)calloc(size1,sizeof(char));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }
   return (tmp);
}

char *alloc_strdup(const char *str) {

   char *tmp=strdup(str);
   if (!tmp)   { printf("strdup failed!\n"); exit(1); }
   return (tmp);
}

double *alloc_double(int size1) {

   double *tmp=(double *)calloc(size1,sizeof(double));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }
   return (tmp);
}

Samples *alloc_sample(int numSample) {

   Samples *tmp=(Samples *)malloc(numSample*sizeof(Samples));
   if (!tmp) { printf("malloc samples failed!\n"); exit(1); }
   return (tmp);
}

Data *alloc_data(int numVariable) {

   Data *tmp=(Data *)malloc(numVariable*sizeof(Data));
   if (!tmp) { printf("malloc data failed!\n"); exit(1); }
   return (tmp);
}


