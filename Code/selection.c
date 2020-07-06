
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "ga_knn.h"


// this subroutine assigns probabilities that are proportional to fitness scores
// Note that this minimizes rather than maximizes
// The minimum is assigned the largest probability
int roulette_wheel_fitness(Fitness *fitness,int populationSize,Wheel *wheel,double convergence) {

   register int i;
   float worst,best,range;

   // fitness sorted in increasing order
   worst=fitness[populationSize-1].value;
   best=fitness[0].value;
   range=worst-best;

   if (range<=convergence) {
      // printf("GA converged ...\n");
      for (i=0; i<populationSize; i++) {
         wheel[i].index=fitness[i].index;
         wheel[i].start=i;
         wheel[i].end=i+1;
      }
      return (1);
   }
   else {
      float scaledScore[populationSize];
      float totalScore=0;
      for (i=0; i<populationSize; i++) {
         scaledScore[i]=worst-fitness[i].value; 
         totalScore += scaledScore[i];
      }
      for (i=0; i<populationSize; i++) scaledScore[i] /= totalScore;

      wheel[0].start=0;
      wheel[0].end  =(float)populationSize*scaledScore[0];
      wheel[0].index=fitness[0].index;

      for (i=1; i<populationSize; i++) {
         wheel[i].start=wheel[i-1].end;
         wheel[i].end=(float)populationSize*scaledScore[i]+wheel[i].start;
         wheel[i].index=fitness[i].index;
      }
      //for (i=0; i<populationSize; i++) printf("%8.5f %8.5f %d\n",wheel[i].start,wheel[i].end,wheel[i].index);
   }
   return (0);
}

