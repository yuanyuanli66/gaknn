#include <stdlib.h>
#include <malloc.h>
#include "ga_knn.h"

void initialize_population(int **allChr,int populationSize,int numVariable,int chromosomeLen) {

   register int m;

   for (m=0; m<populationSize; m++) initialize_chromosome(allChr[m],numVariable,chromosomeLen);
}  

// sample without replacement
void initialize_chromosome(int *chr,int numVariable,int chromosomeLen) {
   
   register int i;
   int numFilled,dummy,found;
   int which[chromosomeLen];
   
   for (i=0; i<chromosomeLen; i++) which[i]=-1;
   numFilled=0;
   while (numFilled<chromosomeLen) {
      dummy=(float)numVariable*genrand_64bits();
      if (dummy==numVariable) dummy--;
      found=0;
      for (i=0; i<numFilled; i++) {
         if (dummy==which[i]) {
            found=1; break; 
         }
      }
      if (!found) {
         which[numFilled]=dummy;
         numFilled++;
      }
   }
   for (i=0; i<chromosomeLen; i++) chr[i]=which[i];
}

