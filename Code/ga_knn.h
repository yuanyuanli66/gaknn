#include <stdio.h>
#include <stdlib.h>

#define  min(a,b)	(((a)<(b))?(a):(b))
#define  max(a,b)	(((a)>(b))?(a):(b))
#define  mod(a,b)	((a)-(b)*((a)/(b)))

#define MAX_CLASS	50
#define MAX_CYCLE	1000
#define LARGE_DISTANCE	999999

typedef struct data_info {
   float *value;
   char *symbol;
   int count;
} Data;

typedef struct sample_info {
   float value;
   char *name;
} Samples;

typedef struct cumulative_prediction {
   float trainingPred[MAX_CYCLE],testingPred[MAX_CYCLE];
   float trainingSum,testingSum;
   int numInTraining,numInTesting;
} CumulativePred;

typedef struct similarity {
   float s;   /* similarity */
   int   p;   /* party      */
} Similarity;

typedef struct score {
   float value;
   int index;
} Fitness;

typedef struct roulette {
   float start,end;
   int index;
} Wheel;

Similarity *alloc_sim(int num);
Wheel *alloc_wheel (int );
Samples *alloc_sample(int );
Data *alloc_data(int );
Fitness *alloc_fitness(int );
CumulativePred *alloc_cum (int size1);

char *alloc_char(int );
char *alloc_strdup(const char *str);
float **alloc_float_float(int ,int );
float *alloc_float(int );
double genrand_64bits(void);

void sort_int(int *data,int size);
void sort_float(float *data,int size);
void sort_fitness(Fitness *,int );
void sort_index(Fitness *fitness,int size);
void sgenrand_64bits(unsigned long long seed);

int *alloc_int (int );
int **alloc_int_int (int ,int );
int read_class(Samples **,const char *);
int read_value(Data **,int ,const char *,Samples *);
void select_training_testing_random(Samples *,int ,int *,int *,float );
int which_chromosome(Wheel *,int );
int num_locus_mutation(int );
int Compare_int(const void *s1, const void *s2);
int Compare_float(const void *s1, const void *s2);
int crossover(int **allChr,int chromosomeLen,int populationSize,Wheel *wheel,Fitness *fitness);
int roulette_wheel_fitness(Fitness *,int ,Wheel *,double );

float calculate_fitness(Samples *,int *,float *,int );

void mutation(int **,int ,int ,int ,Wheel *);
void which_locus(int *,int ,int );
void replace_locus(int ,int *,int ,int *,int );
void initialize_population(int **,int ,int ,int chromosomeLen);
void initialize_chromosome(int *,int ,int );
void predict_training(Samples *,int numTraining,int *,Similarity simTraining[numTraining*numTraining],int ,float *);
void predict_testing(Samples *,int numTraining,int numTesting,int *,Similarity simTesting[numTesting*numTraining],int ,float *);
void print_result(Data *,int ,int *,int ,char *);
void distance_training(const float data_training[/* numVariable*numTraining */],
     int *chr, int chromosomeLen, int numTraining, int *trainingID,
     Similarity simTraining[numTraining*numTraining], int knn);
void distance_training_heap(const float data_training[/* numVariable*numTraining */],
     int *chr, int chromosomeLen, int numTraining, int *trainingID,
     Similarity simTraining[numTraining*numTraining], int knn, float *sum);
void distance_testing(Data *data,int *chr,int chromosomeLen,int numTraining,int *trainingID,int numTesting,
                      int *testingID,Similarity simTesting[numTesting*numTraining],int knn);
void cumulative_prediction(int ,int *,int ,int *,float *,CumulativePred *,Samples *,int );
void print_cumulative(Samples *,int ,CumulativePred *,int ,FILE *,FILE *);
void print_data(Data *data,int ,int );

float spearman_correlation(float *,float *,int );
float pearson_correlation(float *,float *,int );
