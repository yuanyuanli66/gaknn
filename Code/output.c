#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ga_knn.h"

void print_cumulative(Samples *sample,int numSample,CumulativePred *cumulative,int numCycle,FILE *f4,FILE *f5) {

   int i;
   fprintf(f4,"Cycle: %4d\n",numCycle); 
   fprintf(f4,"Sample\tObs.\tPredicted(training)\tPredicted(testing)\n"); 

   for (i=0; i<numSample; i++) {
      if (cumulative[i].numInTraining==0 && cumulative[i].numInTesting==0) {
         fprintf(f4,"%s\t%6.4f\tNA\tNA\tNA\n",sample[i].name,sample[i].value);
      }
      else if (cumulative[i].numInTraining==0 && cumulative[i].numInTesting!=0) {
         fprintf(f4,"%s\t%6.4f\tNA\t%6.4f\n",sample[i].name,sample[i].value,
            cumulative[i].testingSum/(float)cumulative[i].numInTesting);
      }
      else if (cumulative[i].numInTraining!=0 && cumulative[i].numInTesting==0) {
         fprintf(f4,"%s\t%6.4f\t%6.4f\tNA\n",sample[i].name,sample[i].value,
            cumulative[i].trainingSum/(float)cumulative[i].numInTraining);
      }
      else {
         fprintf(f4,"%s\t%6.4f\t",sample[i].name,sample[i].value);
         fprintf(f4,"%6.4f\t%6.4f\n", 
            cumulative[i].trainingSum/(float)cumulative[i].numInTraining,
            cumulative[i].testingSum/(float)cumulative[i].numInTesting);
      }
   }
   fflush(f4); 

   float predTr[numSample],obsdTr[numSample];
   float predTe[numSample],obsdTe[numSample];
   int cnTr=0,cnTe=0;
   for (i=0; i<numSample; i++) {
      if (cumulative[i].numInTraining!=0 && sample[i].value !=-9999) {
         predTr[cnTr]=cumulative[i].trainingSum/(float)cumulative[i].numInTraining;
         obsdTr[cnTr]=sample[i].value; 
         cnTr++; 
      }
      if (cumulative[i].numInTesting!=0 && sample[i].value !=-9999) {
         predTe[cnTe]=cumulative[i].testingSum/(float)cumulative[i].numInTesting;
         obsdTe[cnTe]=sample[i].value; 
         cnTe++; 
      }
   }
   float rmse_tr=9999,rmse_te=9999;
   float p_tr=-2,s_tr=-2,p_te=-2,s_te=-2;
   if (cnTr>1) {
      rmse_tr=0; for (i=0; i<cnTr; i++) rmse_tr += (predTr[i]-obsdTr[i])*(predTr[i]-obsdTr[i]); rmse_tr=sqrt(rmse_tr/(float)(cnTr-1));
      p_tr=pearson_correlation(predTr,obsdTr,cnTr);
      s_tr=spearman_correlation(predTr,obsdTr,cnTr);
   }
   if (cnTe>1) {
      rmse_te=0; for (i=0; i<cnTe; i++) rmse_te += (predTe[i]-obsdTe[i])*(predTe[i]-obsdTe[i]); rmse_te=sqrt(rmse_te/(float)(cnTe-1));
      p_te=pearson_correlation(predTe,obsdTe,cnTe);
      s_te=spearman_correlation(predTe,obsdTe,cnTe);
   }
   fprintf(f5,"Cycle: %d\t",numCycle);
   fprintf(f5,"%7.5f\t%5.4f\t%5.4f\t%7.5f\t%5.4f\t%5.4f\n",rmse_tr,p_tr,s_tr,rmse_te,p_te,s_te);
   fflush(f5); 
}


void cumulative_prediction(int numTraining,int *trainingID,int numTesting,int *testingID,
   float *predictedValue,CumulativePred *cumulative,Samples *sample,int cycle) {

   register int i;
   int wh;

   for (i=0; i<numTraining; i++) {
      wh=trainingID[i];
      cumulative[wh].trainingSum +=predictedValue[wh];
      cumulative[wh].trainingPred[cycle]=predictedValue[wh];
      (cumulative[wh].numInTraining)++;
   }

   for (i=0; i<numTesting; i++) {
      wh=testingID[i];
      cumulative[wh].testingSum +=predictedValue[wh];
      cumulative[wh].testingPred[cycle]=predictedValue[wh];
      (cumulative[wh].numInTesting)++;
   }
}

void print_result(Data *data,int numVariable,int *trainingID,int numTraining,char *outFileCount) {

   FILE *f2;
   register int i;

   f2=fopen(outFileCount,"w");
   for (i=0; i<numVariable; i++) {
      fprintf(f2,"%s\t%d\n",data[i].symbol,data[i].count);
   }
   fclose(f2);
}

void print_data(Data *data,int numVariable,int numSample) {

   FILE *fp;
   int i;

   fp=fopen("debug.data","w");

   fprintf(fp,"%s\t",data[0].symbol);
   for (i=0; i<numSample; i++) fprintf(fp,"%5.4f\t",data[0].value[i]); fprintf(fp,"\n");

   fprintf(fp,"%s\t",data[numVariable-1].symbol);
   for (i=0; i<numSample; i++) fprintf(fp,"%5.4f\t",data[numVariable-1].value[i]); fprintf(fp,"\n");

   fclose(fp);
}

