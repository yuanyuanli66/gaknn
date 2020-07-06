
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "ga_knn.h"
#define	GIVE_UP	0

void sort_int(int *,int );

int crossover(int **chromosome,int chromosomeLen,int populationSize,Wheel *wheel,Fitness *fitness) {

   int wh1,wh2,crossPoint,du;
   int inChr1,inChr2;
   int count;
   int newChr[populationSize][chromosomeLen],tmpChr1[chromosomeLen],tmpChr2[chromosomeLen];
   register int i,j;

   int numSaved=1;
   for (i=0; i<numSaved; i++) {
      for (j=0; j<chromosomeLen; j++) newChr[i][j]=chromosome[wheel[i].index][j];
   }
   count=numSaved;
   do {
      // Pick two chromosomes
      int numTry=0;
      do {
         wh1=which_chromosome(wheel,populationSize);
         wh2=which_chromosome(wheel,populationSize);
         numTry++;
         if (numTry==100) return (GIVE_UP); 
      } while (fitness[wh1].value == fitness[wh2].value);

      for (i=0; i<chromosomeLen; i++) tmpChr1[i]=chromosome[wh1][i];
      for (i=0; i<chromosomeLen; i++) tmpChr2[i]=chromosome[wh2][i];

      sort_int(tmpChr1,chromosomeLen); // increasing order for easy checking overlaps
      sort_int(tmpChr2,chromosomeLen); // same as above

      // find a crossover point
      crossPoint=(int)((float)chromosomeLen*genrand_64bits()); 
      if (crossPoint==0) crossPoint++;
      if (crossPoint==chromosomeLen-1) crossPoint--;
   
      // exchange elements (up to crossover point) between two chromosomes 
      for (i=0; i<crossPoint; i++) {
         du=tmpChr1[i];
         tmpChr1[i]=tmpChr2[i];
         tmpChr2[i]=du; 
      }
      // check for duplicate elements in a chromosome
      inChr1=0; inChr2=0;
      for (i=0; i<crossPoint; i++) {
         for (j=crossPoint; j<chromosomeLen; j++) {
            if (tmpChr1[i]==tmpChr1[j]) { inChr1=1; break; } 
         } 
         for (j=crossPoint; j<chromosomeLen; j++) {
            if (tmpChr2[i]==tmpChr2[j]) { inChr2=1; break; } 
         } 
      }
      // if no overlaps, add to the population and update count
      if (!inChr1) {
         for (i=0; i<chromosomeLen; i++) newChr[count][i]=tmpChr1[i]; 
         count++;
         if (count==populationSize) break;  
      }   
      if (!inChr2) {
         for (i=0; i<chromosomeLen; i++) newChr[count][i]=tmpChr2[i]; 
         count++;
         if (count==populationSize) break;  
      }   
   } while (count!=populationSize);

   for (i=0; i<populationSize; i++) {
      for (j=0; j<chromosomeLen; j++) chromosome[i][j]=newChr[i][j];
   }
   return 1;
}

