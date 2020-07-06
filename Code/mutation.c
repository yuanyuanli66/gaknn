
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "ga_knn.h"

void mutation(int **chromosome,int numVariable,int chromosomeLen,int populationSize,Wheel *wheel) {

   int wh,numLocusMutate;
   int newChr[populationSize][chromosomeLen],tmpChr[chromosomeLen];
   int locus[chromosomeLen];
   register int i,j;

   // save the best
   int numSaved=1;
   for (i=0; i<numSaved; i++) {
      for (j=0; j<chromosomeLen; j++) newChr[i][j]=chromosome[wheel[i].index][j];
   }

   int count=numSaved; 
   do {
      // pick a chromosome
      wh=which_chromosome(wheel,populationSize);
      for (i=0; i<chromosomeLen; i++) tmpChr[i]=chromosome[wh][i];

      // number of loci to be mutated 
      numLocusMutate=num_locus_mutation(chromosomeLen);

      // pick loci
      which_locus(locus,chromosomeLen,numLocusMutate);
    
      // replace loci 
      replace_locus(numVariable,tmpChr,chromosomeLen,locus,numLocusMutate);
      for (i=0; i<chromosomeLen; i++) newChr[count][i]=tmpChr[i]; 

      count++;
   } while (count<(int)(0.99*populationSize));

   for (i=count; i<populationSize; i++) initialize_chromosome(newChr[i],numVariable,chromosomeLen);

   for (i=0; i<populationSize; i++) {
      for (j=0; j<chromosomeLen; j++) chromosome[i][j]=newChr[i][j];
   } 
}

int which_chromosome(Wheel *wheel,int populationSize) {

   register int i;
   int wh;
   float dummy;

   dummy=(float)populationSize*genrand_64bits();
   for (i=0; i<populationSize; i++) {
      if (dummy>=wheel[i].start && dummy<=wheel[i].end) return (wheel[i].index); 
   }
   wh=(int)((float)populationSize*genrand_64bits());
   if (wh==populationSize) wh--;
   return (wh);
}

int num_locus_mutation(int chromosomeLen) {
 
   float r;
   register int i;

   r=genrand_64bits();
   for (i=1; i<=chromosomeLen; i++) {
      if (r>=1.0/pow(2.0,(float)i))  return (i);
   }
   return (1);
}

void which_locus(int *locus,int chromosomeLen,int numLocusMutate) {

   register int i;
   int numFilled,dummy,found;

   for (i=0; i<numLocusMutate; i++) locus[i]=-1;

   numFilled=0;
   while (numFilled<numLocusMutate) {
      dummy=(int)((float)chromosomeLen*genrand_64bits());
      if (dummy==chromosomeLen) dummy--;

      found=0;
      for (i=0; i<numFilled; i++) {
         if (dummy==locus[i]) {
            found=1; break;
         }
      }
      if (!found) {
         locus[numFilled]=dummy;
         numFilled++;
      }
   }
}

void replace_locus(int numVariable,int *tmpChr,int chromosomeLen,int *locus,int numLocusMutate) {

   register int j;
   int dummy;

   int cn=0;
   do {
      dummy=(int)((float)numVariable*genrand_64bits());
      if (dummy==numVariable) dummy--;
      int inChr=0;
      for (j=0; j<chromosomeLen; j++) {
         if (dummy==tmpChr[j]) {
            inChr=1; break;
         }
      }
      if (!inChr) { tmpChr[locus[cn]]=dummy; cn++; }
   } while (cn<numLocusMutate);
}

