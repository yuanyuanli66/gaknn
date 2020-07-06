/*--------------------------------------------------------------
* Sample classification and gene assessment for expression data *
*   using a genetic algorithm and k-nearest neighbor method     *
*                                                               *
*                        Leping Li                              *
*        National Institute of Environmental Health Sciences    *
*                  National Institute of Health                 *
*                                                               *
*                  Date: March 26, 1999                         *
*     Last Modification: July  19, 2001                         *
*              Revision: December 25, 2013                      *
*         last revision: July 04, 2020                          *
----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ga_knn.h"

int main(int argc,char **argv) {

   Samples *sample;                // sample name and associated value, e.g., IC50 if known, otherwised labeled as -9999
   Data *data;                     // expression data - all
   float *data_training;           // expression data - training 
   Wheel *wheel;	   	   // roulette wheel selection of chromosomes in population
   Similarity *simTesting;         // 2D array (numTesting,numTraining)
   CumulativePred *cumulative;	   // all samples cumulative prediction results

   int numSample;
   int numTraining,numTesting; 
   int numVariable;                // number of genes
   int **allChr;                   // store all candidate chromosomes in GA population
   int *trainingID,*testingID;
   float  **predicted;             // predicted values for each chromosome - matrix
   float *predictedValue;          // predicted values for the "best" chromosome - vector length of sample
   float *score;                   // temporary score to keep track of convergence
   Fitness *fitness;               // goes with wheel for chromosomes
   time_t start,finish;
   clock_t cpu_start, cpu_finish;
   struct timespec ts_start, ts_finish;
   double cpu_time, wall_time;
   register int i,j,cycle,gen;
   FILE *fp,*f2,*f4,*f5;

   if (argc==1) {
      printf("\n\nusage: %s -classFile classFileName -dataFile dataFileName\n"
             "       \nadditional\n"
             "       -knn\t\tinteger[3, 5,or 7] (default 5)\n"
             "       -popSize\t\tinteger (default 2000)\n"
             "       -chromosome\tinteger(default 30)\n"
             "       -numGen\t\tinteger (default 1000)\n"
             "       -numCycle\tinteger (default 100)\n"
             "       -propTest\t[0 to 0.5] (default 0.20)\n"
             "       -thread\t\t[>=1,-1=nproc] (default 50)\n"
             "       -seed\t\tinteger (default time(0)\n"
             "       -step\t\tinteger(default 300)\n"
             "       -mutProb\t\t[0-1](default 0.3333)\n" 
             "       -nprint\t\tinteger(default 100)\n\n",
             argv[0]);
      exit(EXIT_SUCCESS);
   }

   // defaults:
   int knn=5;
   int chromosomeLen=30;
   int populationSize=2000;
   int numGen=1000;
   int numCycle=100;
   int numThreads=50;
   float propTest=0.20;
   int use_heap=0;
   int nprint=100;
   int stepNoChange=300;             // if no change in fitness value, 'converged'
   float mutProb=0.3333333;
   double convergence=0.001;
   const char *classFileName=NULL;
   const char *dataFileName=NULL;
   char *outFileInfo=NULL; 
   char *outFileCount=NULL;
   char *outFileChr=NULL; 
   char *outFilePred=NULL; 
   char *outFileAccuracy=NULL; 
   unsigned long long seed=time(0);

   // default output file names in the directory in which the program is running
   outFileChr=strdup("top_chromosome.txt");
   outFileCount=strdup("selection_count.txt");
   outFileInfo=strdup("info.txt");
   outFilePred=strdup("cumulative_prediction.txt");
   outFileAccuracy=strdup("accuracy.txt");

   for (i=1; i<argc; i++) {
      if      (!strcmp(argv[i],"-classFile") && i<argc-1)   classFileName=argv[++i]; 
      else if (!strcmp(argv[i],"-dataFile") && i<argc-1)    dataFileName=argv[++i]; 
      else if (!strcmp(argv[i],"-knn") && i<argc-1)         knn=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-popSize") && i<argc-1)     populationSize=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-chromosome") && i<argc-1)  chromosomeLen=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-numGen") && i<argc-1)      numGen=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-numCycle") && i<argc-1)    numCycle=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-propTest") && i<argc-1)    propTest=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-convergence") && i<argc-1) convergence=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-mutProb") && i<argc-1)     mutProb=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-thread") && i<argc-1)      numThreads=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-seed") && i<argc-1)        seed=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-nprint") && i<argc-1)      nprint=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-step") && i<argc-1)        stepNoChange=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-outInfo") && i<argc-1)     outFileInfo=strdup(argv[++i]); 
      else if (!strcmp(argv[i],"-outChr") && i<argc-1)      outFileChr=strdup(argv[++i]); 
      else if (!strcmp(argv[i],"-outCount") && i<argc-1)    outFileCount=strdup(argv[++i]); 
      else if (!strcmp(argv[i],"-outPred") && i<argc-1)     outFilePred=strdup(argv[++i]); 
      else if (!strcmp(argv[i],"-outAccuracy") && i<argc-1) outFileAccuracy=strdup(argv[++i]); 
      else if (!strcmp(argv[i],"-heap"))                    use_heap=1; 
      else { printf("argument: %s unknown\n",argv[i]); exit(1); }
   }

   if (populationSize<2) { printf("Population size [%d] too small\n",populationSize); exit(1); }
   if (propTest==0) printf("all samples are used in training...\n\n");

   if (numThreads==0) { printf("number of Threads==0, reset to 1\n"); numThreads=1; }
#ifndef _OPENMP
   if (numThreads!=1) { printf("Non-parallel code, number of threads ignored\n"); numThreads=1; }
#endif

   if (!classFileName) { printf("Error: -classFile not defined\n"); exit(1); }
   numSample=read_class(&sample,classFileName);
   // FILE *fd=fopen("debug.class","w");
   //for (i=0; i<numSample; i++) fprintf(fd,"%s\t%5.5f\n",sample[i].name,sample[i].value); 
   //fclose(fd);

   if (!dataFileName) { printf("Error: -dataFile not defined\n"); exit(1); }
   numVariable=read_value(&data,numSample,dataFileName,sample);

   printf("number of samples:\t\t%d\n",numSample);
   printf("number of variables:\t\t%d\n",numVariable);
   printf("percentage of testing:\t\t%4.2f\n",propTest);
   printf("population size:\t\t%d\n",populationSize);
   printf("number of GA generations:\t%d\n",numGen);
   printf("number of cycles set:\t\t%d\n",numCycle);
   printf("number of threads:\t\t%d\n",numThreads);

   sgenrand_64bits(seed); // before training and testing split

   trainingID=alloc_int(numSample);
   testingID=alloc_int(numSample);

   //---------------------------------------------------
   // find out which samples with known outcome values. 
   // those to be predicted must have values of -9999 
   //---------------------------------------------------
   int unknownSampleCn=0;
   for (i=0; i<numSample; i++) {
      if (sample[i].value ==-9999) unknownSampleCn++; 
   }
   int knownSampleCn=numSample-unknownSampleCn;
   printf("\nNumber of samples with known IC50 values:\t%d\n",knownSampleCn);

   if (knownSampleCn==numSample) printf("performing Monte Carlo sampling and cross-validation using all samples\n");
   else                          printf("performing Monte Carlo sampling using samples with known IC50s and predict unknown samples\n");

   numTraining=knownSampleCn*(1.0-propTest);
   numTesting=numSample-numTraining;
 
   printf("number training and testing: %d %d\n",numTraining,numTesting);

   // Check stack limit and increase if necessary, for stack-allocated arrays.
   // simTraining[] and sum[]
   if (!use_heap) {
      struct rlimit rlim;
      rlim_t need=numTraining * numTraining * sizeof(struct similarity); // simTraining[]
      need += sizeof(float)*numTraining*(numTraining-1)/2; // temporary sum[] array
      need += 2048 * 1024; // 2M for other stack usage (more than generous)
      printf("estimated min. stack:\t%luK\n",(unsigned long)(need/1024));
#ifdef _OPENMP
      char size[32];
      snprintf(size,sizeof(size),"%luK",(unsigned long)need/1024);
      printf("setenv OMP_STACKSIZE=%s\n",size);
      setenv("OMP_STACKSIZE",size,0);
#endif
      if (getrlimit(RLIMIT_STACK, &rlim)) {
        printf("Warning: Cannot read stacksize limit.\n");
      } else {
         if (rlim.rlim_cur != RLIM_INFINITY && rlim.rlim_cur < need) {
           printf("Increasing stack limit to %luK\n",(unsigned long)(need/1024));
           rlim.rlim_cur=need;
           if (setrlimit(RLIMIT_STACK,&rlim)) {
              printf("Warning: Unable to increase stack size limit, using heap allocation.\n");
              use_heap=1;
           }
         } else {
           printf("execution stack size:\t");
           if (rlim.rlim_cur == RLIM_INFINITY) puts("infinity");
           else printf("%luK\n",(unsigned long)(rlim.rlim_cur/1024));
         }
      }
   }
#ifdef _OPENMP
   if (numThreads<0) numThreads=omp_get_max_threads();
   else omp_set_num_threads(numThreads);
#endif

   allChr=alloc_int_int(populationSize,chromosomeLen);
   wheel=alloc_wheel(populationSize);
   fitness=alloc_fitness(populationSize);

   simTesting=alloc_sim(numTesting*numTraining);
   cumulative=alloc_cum(numSample);
   predictedValue=alloc_float(numSample);
   predicted=alloc_float_float(populationSize,numSample);
   score=alloc_float(numGen);

   fp=fopen(outFileInfo,"w"); 
   f2=fopen(outFileChr,"w");
   
   time(&start);
   fprintf(fp,"##################################################################\n");
   fprintf(fp,"         Last update: May 30, 2019\n\n");
   fprintf(fp,"Random number generator: Marsaglia's algorithm.\n");
   fprintf(fp,"Random seed:\t\t\t\t\t%lld\n",seed);
   fprintf(fp,"Processor ID:\t\t\t\t\t%d\n",getpid());
   fprintf(fp,"Number of threads specified:\t\t\t%d\n",numThreads);
   fprintf(fp,"Job started: %s\n",asctime(localtime(&start)));

   fprintf(fp,"Parameters:\n");
   fprintf(fp,"  Number of cycles:\t\t\t\t%d\n",numCycle);
   fprintf(fp,"  Number of generations & population size:\t%d, %d\n",numGen,populationSize);
   fprintf(fp,"  Chromosome length (feature dimension):\t%d\n",chromosomeLen);
   fprintf(fp,"  Mutation probability (crossover=1-that):\t%4.3f\n",mutProb);
   fprintf(fp,"  Convergence:\n");
   fprintf(fp,"     fitness score change<=%6.5f for %d generations\n",convergence,stepNoChange);
   fprintf(fp,"  KNN:\n");
   fprintf(fp,"     K-nearest neighbors:\t\t\t%2d\n",knn);
   fprintf(fp,"\nData info:\n");
   fprintf(fp,"  Number of samples:\t\t\t\t%d\n",numSample);
   if (knownSampleCn==numSample) fprintf(fp, "      Performing Monte Carlo sampling and cross-validation using all samples\n");
   else                          fprintf(fp, "      Performing Monte Carlo sampling using samples with known outcome values\n");
   fprintf(fp,"  Numbers samples in training and testing:\t%d, %d\n",numTraining,numTesting);
   fprintf(fp,"  Number of variables:\t\t\t\t%d\n",numVariable);
   fprintf(fp,"  Data file:\t%s\n",dataFileName);
   fprintf(fp,"  Calss file:\t%s\n",classFileName);
   fprintf(fp,"##################################################################\n");
   fflush (fp);

   data_training=alloc_float(numVariable * numTraining);

   for (i=0; i<numSample; i++) {
      cumulative[i].trainingSum=0;
      cumulative[i].testingSum=0;
      cumulative[i].numInTraining=0;
      cumulative[i].numInTesting=0;
   }
   for (i=0; i<numVariable; i++) data[i].count=0;

   cpu_start=clock();
   // Could use omp_get_wtime() here, with OpenMP.
   clock_gettime(CLOCK_REALTIME,&ts_start);
   // repeat the GA 'numCycle' times

   f5=fopen(outFileAccuracy,"w");

   for (cycle=0; cycle<numCycle; cycle++) {

      which_locus(trainingID,knownSampleCn,numTraining);
      int cn=0;
      for (j=0; j<numSample; j++) {
         int found=0;
         for (i=0; i<numTraining; i++) {
            if (j==trainingID[i]) { found=1; break; } 
         }
         if (!found) {
            testingID[cn]=j; cn++; 
         }
      }
 
      // select_training_testing_random(sample,numSample,trainingID,testingID,propTest);
      for (i=0; i<numVariable; i++) {
         for (j=0; j<numTraining; j++) data_training[i*numTraining+j]=data[i].value[trainingID[j]];
      }

      initialize_population(allChr,populationSize,numVariable,chromosomeLen);
      printf("Done - initializing chromosomes. cycle=%3d/%d;",cycle+1,numCycle);

      // a single GA run
      for (gen=0; gen<numGen; gen++) {

      // process each "chromosome" in the population 
      // 1) compute the Euclidean distance between a pair of samples in the chromosomes dimension
      // 2) find the k-nearest neighbors
      // 3) average the k training outcomes as the predicted value for the sample
      // 4) compute the fitness score for each chromosome in the population

      /********* Stack version *********/
         if (!use_heap)
#pragma omp parallel
         {
            int id;
            // Automatic array, allocated on the stack, gives better performance than heap arrays.
            Similarity simTraining[numTraining*numTraining];
            // temporary sum[] array is allocated on the stack in distance_training()
#pragma omp  for
             for (id=0; id<populationSize; id++) {
               distance_training(data_training,allChr[id],chromosomeLen,numTraining,trainingID,simTraining,knn);
               predict_training(sample,numTraining,trainingID,simTraining,knn,predicted[id]);
               fitness[id].value=calculate_fitness(sample,trainingID,predicted[id],numTraining);
               fitness[id].index=id;
            }
         }
/********* Heap version *********/
         else
#pragma omp parallel
         {
            int id;
            Similarity *simTraining=alloc_sim(numTraining*numTraining);
            float *sum=malloc(sizeof(float)*numTraining*(numTraining-1)/2);
            if (!sum) { fprintf(stderr,"Error allocating working sum\n"); exit(1); }
#pragma omp  for
             for (id=0; id<populationSize; id++) {
               distance_training_heap(data_training,allChr[id],chromosomeLen,numTraining,trainingID,simTraining,knn,sum);
               predict_training(sample,numTraining,trainingID,simTraining,knn,predicted[id]);
               fitness[id].value=calculate_fitness(sample,trainingID,predicted[id],numTraining);
               fitness[id].index=id;
            }
            free(sum);
            free(simTraining);
         }
/*********************************/

         sort_fitness(fitness,populationSize);
         //printf("cycle[%4d] generation[%4d] fitness: %5.3f\n",cycle+1,gen+1,fitness[0].value);
         if (gen%nprint==0 || gen==numGen-1) {
            fprintf(fp,"cycle[%4d] generation[%4d] fitness: %5.3f\n",cycle+1,gen+1,fitness[0].value);
            fflush(fp);
         }

         for (i=0; i<numTraining; i++) predictedValue[trainingID[i]]=predicted[fitness[0].index][trainingID[i]];

         score[gen]=fitness[0].value;
         if (gen+1>min(stepNoChange+1,numGen)) {
            if (fabsf(score[gen]-score[gen-stepNoChange-1])<convergence) { 
               printf("converged...\n");
               fprintf(fp,"cycle[%4d] generation[%4d] fitness: %5.3f\n",cycle+1,gen+1,fitness[0].value);
               fflush(fp);
               break;
            }
         }

         // roulette wheel selection of new generation based on the fitness of the previous generation
         //  -- survival of the fittest
         int homogeneous=roulette_wheel_fitness(fitness,populationSize,wheel,convergence);
         // "genetic" mutation of the chromosomes - replacing a gene(s) on the chromosome
         if (homogeneous) {
            // printf("population homogeneous... break...\n");
            mutation(allChr,numVariable,chromosomeLen,populationSize,wheel);
         }
         else {
            if (genrand_64bits()<=mutProb) { mutation(allChr,numVariable,chromosomeLen,populationSize,wheel); }
            else { 
               int success=1;
               success=crossover(allChr,chromosomeLen,populationSize,wheel,fitness); 
               //printf("crossover - %d.\n",success);
               if (!success) mutation(allChr,numVariable,chromosomeLen,populationSize,wheel); 
            }
         }
      } // end gen

      for (i=0; i<chromosomeLen; i++) data[allChr[0][i]].count++;
      print_result(data,numVariable,trainingID,numTraining,outFileCount);
   
      // compute the Euclidean distance between each test sample and each training sample 
      distance_testing(data,allChr[0],chromosomeLen,numTraining,trainingID,numTesting,testingID,simTesting,knn);
      // predict test samples based on the classes of their nearest neighbors of the training samples 
      predict_testing(sample,numTraining,numTesting,testingID,simTesting,knn,predictedValue);
      // summarize the results from multiple GA runs
      cumulative_prediction(numTraining,trainingID,numTesting,testingID,predictedValue,cumulative,sample,cycle); 

      f4=fopen(outFilePred,"w");
      print_cumulative(sample,numSample,cumulative,cycle+1,f4,f5);
      fclose(f4);

      for (i=0; i<chromosomeLen; i++) fprintf(f2,"%d\t",allChr[0][i]); 
      fprintf(f2,"\t%6.3f\n",fitness[0].value);
      fflush(f2);

      clock_gettime(CLOCK_REALTIME,&ts_finish);
      wall_time=((double)(ts_finish.tv_sec)+(double)(ts_finish.tv_nsec)/1e+9) -
           ((double)(ts_start.tv_sec)+(double)(ts_start.tv_nsec)/1e+9);
      printf("\t%.3f sec/cycle\n",wall_time/(double)(cycle+1));

   } // end cycle
   fclose(f2); fclose(f5);
   time(&finish);
   fprintf(fp,"job finished:\t%s\n",asctime(localtime(&finish)));
   cpu_finish=clock();
   cpu_time=((double) (cpu_finish - cpu_start)) / CLOCKS_PER_SEC;
   clock_gettime(CLOCK_REALTIME,&ts_finish);
   wall_time=((double)(ts_finish.tv_sec)+(double)(ts_finish.tv_nsec)/1e+9) -
        ((double)(ts_start.tv_sec)+(double)(ts_start.tv_nsec)/1e+9);
   fprintf(fp,"CPU time: %.2f sec\n",cpu_time);
   fprintf(fp,"Wall time: %.2f sec\n",wall_time);
   fprintf(fp,"Efficiency: %.2f%%\n",100.0*cpu_time/(double)numThreads/wall_time);
   printf("CPU time: %.2f sec\n",cpu_time);
   printf("Wall time: %.2f sec\n",wall_time);
   printf("Efficiency: %.2f%%\n",100.0*cpu_time/(double)numThreads/wall_time);
   fclose(fp);

   return (EXIT_SUCCESS);
}

