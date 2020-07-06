#include <stdlib.h>
#include "ga_knn.h"

// increasing order
static int Compare_index(const void *s1, const void *s2) {
  return ((Fitness *)s1)->index - ((Fitness *)s2)->index;
}

/* qsort */
void sort_index(Fitness *fitness,int size) {
   qsort((void *)fitness,(size_t)size,sizeof(Fitness),Compare_index);
}

// increasing order
static int Compare_fitness(const void *s1, const void *s2) {
   if (((Fitness *)s1)->value < ((Fitness *)s2)->value) { return -1; }
   if (((Fitness *)s1)->value > ((Fitness *)s2)->value) { return  1; }
   return 0;
}

/* qsort */
void sort_fitness(Fitness *fitness,int size) {
   qsort((void *)fitness,(size_t)size,sizeof(Fitness),Compare_fitness);
}

/* increasing */
void sort_int(int *data,int size) {
   int (*compar)(const void *,const void *);
   compar = Compare_int;
   qsort((void *)data,(size_t)size,sizeof(int),compar);
}

int Compare_int(const void *s1, const void *s2) {
   if (*((int *)s1)< (*(int *)s2)) { return -1; }
   if (*((int *)s1)> (*(int *)s2)) { return  1; }
      return 0;
}

/* increasing */
void sort_float(float *data,int size) {
   int (*compar)(const void *,const void *);
   compar = Compare_float;
   qsort((void *)data,(size_t)size,sizeof(int),compar);
}

int Compare_float(const void *s1, const void *s2) {
   if (*((float *)s1)< (*(float *)s2)) { return -1; }
   if (*((float *)s1)> (*(float *)s2)) { return  1; }
      return 0;
}

