#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "ga_knn.h"

// 2D array macros, to make code slightly more readable.
// Both sim arrays have the inner dimension 'numTraining'
#define simTraining_(i,j) (simTraining[i*numTraining+j])
#define simTesting_(i,j) (simTesting[i*numTraining+j])

void predict_training(Samples *sample,int numTraining,int *trainingID,Similarity simTraining[numTraining*numTraining],int knn,float *predictedValue) {

   int i,j;
   for (i=0; i<numTraining; i++) {
      float ave=0;
      for (j=0; j<knn; j++) {
         int which=simTraining_(i,j).p; // training j nearest training neighbors 
         ave +=sample[which].value;
      }
      predictedValue[trainingID[i]]=ave/(float)knn;
   }
}

void predict_testing(Samples *sample,int numTraining,int numTesting,int *testingID,Similarity simTesting[numTesting*numTraining],int knn,float *predictedValue) {

   int i,j;
   for (i=0; i<numTesting; i++) {
      float ave=0;
      for (j=0; j<knn; j++) {
         int which=simTesting_(i,j).p; // testing j nearest training neighbors 
         ave += sample[which].value;
      }
      predictedValue[testingID[i]]=ave/(float)knn;
   }
}

// Nearly all CPU time is spent in partial_sort_ascen() and distance_training().
// With current optimizations, the sort procedure uses most of the CPU time.

// Insertion sort, with the lower sorted region limited to KNN
// Unsorted data is not maintained. (Some values are simply dicarded.)
static inline void partial_sort_ascen(Similarity x[/* size */],int size,int knn) {
   int i,j;
   // Standard insertion sort of the initial KNN elements
   for (i=1;i<knn;i++) {
      Similarity tmp = x[i];
      for (j=i; j>0 && x[j-1].s > tmp.s; j--) x[j]=x[j-1];
      x[j] = tmp;
   }
   // Truncated insertion sort of the remaining >KNN elements
   // All remaining elements are sorted to the first KNN in a single pass.
   for (i=knn;i<size;i++) {
      // Only shift in if it's lower than the largest of the KNN subgroup. 
      if (x[knn-1].s > x[i].s) {
         // Shift larger values upwards; element KNN is simply overwritten.
         for (j=knn-1; j>0 && x[j-1].s > x[i].s; j--) x[j]=x[j-1];
         // Store new value in the correct sort position.
         x[j] = x[i];
      }
   }
}

// data[] lookup (k loop) moved to the outermost loop, added local temporary sum[] array
void distance_training(const float data_training[/* numVariable*numTraining */],
     int *chr, int chromosomeLen, int numTraining, int *trainingID,
     Similarity simTraining[numTraining*numTraining], int knn) {
   register int i,j,k,n;
   float sum[numTraining*(numTraining-1)/2];

   // Initialize the temporary local sum array
   bzero(sum,sizeof(float)*numTraining*(numTraining-1)/2);

   // compute upper echelon (CPU intensive loop)
   // Temporaries '*p' and 'pi' should not be necessary, (should be automatic
   // by optimizer) but seems to make vectorizer more effective.
   for (k=0; k<chromosomeLen; k++) {
      const float *p = &data_training[chr[k]*numTraining];
      for (i=0,n=0; i<numTraining-1; i++) {
         const float pi = p[i];
         for (j=i+1; j<numTraining; j++,n++) {
            float diff = pi - p[j];
            sum[n] += diff*diff;
         }
      }
   }

    // initialize diagonal and copy sums to array
   for (i=0,n=0; i<numTraining; i++) {
      simTraining_(i,i).s=LARGE_DISTANCE;
      simTraining_(i,i).p=trainingID[i];
      for (j=i+1; j<numTraining; j++,n++) {
         simTraining_(j,i).s=simTraining_(i,j).s=sum[n];
         simTraining_(j,i).p=trainingID[i];
         simTraining_(i,j).p=trainingID[j];
      } 
   }

   // sort results
   for (i=0; i<numTraining; i++) partial_sort_ascen(&simTraining_(i,0),numTraining,knn);
}

// Version with sum[] array heap-allocated
void distance_training_heap(const float data_training[/* numVariable*numTraining */],
     int *chr, int chromosomeLen, int numTraining, int *trainingID,
     Similarity simTraining[numTraining*numTraining], int knn, float *sum) {
   register int i,j,k,n;

   bzero(sum,sizeof(float)*numTraining*(numTraining-1)/2);

   for (k=0; k<chromosomeLen; k++) {
      const float *p = &data_training[chr[k]*numTraining];
      for (i=0,n=0; i<numTraining-1; i++) {
         float pi = p[i];
         for (j=i+1; j<numTraining; j++,n++) {
            float diff = pi - p[j];
            sum[n] += diff*diff;
         }
      }
   }

   for (i=0,n=0; i<numTraining; i++) {
      simTraining_(i,i).s=LARGE_DISTANCE;
      simTraining_(i,i).p=trainingID[i];
      for (j=i+1; j<numTraining; j++,n++) {
         simTraining_(j,i).s=simTraining_(i,j).s=sum[n];
         simTraining_(j,i).p=trainingID[i];
         simTraining_(i,j).p=trainingID[j];
      } 
   }

   for (i=0; i<numTraining; i++) partial_sort_ascen(&simTraining_(i,0),numTraining,knn);
}

void distance_testing(Data *data,int *chr,int chromosomeLen,int numTraining,int *trainingID,int numTesting,
                      int *testingID,Similarity simTesting[numTesting*numTraining],int knn) {
#pragma omp parallel
   {
      int i;
#pragma omp for
      for (i=0; i<numTesting; i++) {
         int j,k;
         for (j=0; j<numTraining; j++) {
            float sum=0.0; 
            for (k=0; k<chromosomeLen; k++) {
               float diff=data[chr[k]].value[testingID[i]]-data[chr[k]].value[trainingID[j]];
               sum += diff*diff;
            }
            simTesting_(i,j).s=sum;
            simTesting_(i,j).p=trainingID[j];
         }
         partial_sort_ascen(&simTesting_(i,0),numTraining,knn+1);
      }
   }
}

float calculate_fitness(Samples *sample,int *trainingID,float *predictedValue,int numTraining) {

   register int i;
   float diff,sum;

   sum=0;
   for (i=0; i<numTraining; i++) {
      diff=predictedValue[trainingID[i]]-sample[trainingID[i]].value;
      sum +=diff*diff;
   }
   return (sum);
}

