#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ga_knn.h"

// Read samples into a temporary linked-list, then transfer into a flat array.
int read_class(Samples **sample,const char *fileName) {

   FILE *fp;
   int i,cn,len;
   char buffer[5000],*tok;
   struct sample_info_list {
      float value;
      char *name;
      struct sample_info_list *next;
   } *p_first=NULL, *p_last, *p_sample;

   fp=fopen(fileName,"r");
   if (!fp) { perror(fileName); exit(1); }

   printf("reading... %s\n",fileName);

   int foundUNK=0;
   for (cn=0; !feof(fp); cn++) {
      int unk=0;
      if (!fgets(buffer,sizeof(buffer),fp)) break;
      len=strlen(buffer);
      buffer[len-1]='\0';
      p_sample=malloc(sizeof(*p_sample));
      if (!p_sample) { printf("Error allocating sample!\n"); exit(1); }
      if (!p_first) p_first=p_sample;
      else p_last->next=p_sample;
      p_last=p_sample;
      tok=strtok(buffer,"\t");
      p_sample->name=alloc_strdup(tok);
      tok=strtok(0,"\t");
      if (tok[0]=='N' && tok[1]=='A') p_sample->value=-9999;
      else                            p_sample->value=atof(tok);
      if (p_sample->value==-9999) {
         unk=1; foundUNK=1; 
      }
      if (foundUNK==1 && unk==0) { 
         printf("\nError: the unknown samples must be placed after the known samples\n"); 
         printf("Example:\n");
         printf("sample1\t0.58\n");
         printf("sample2\t1.22\n");
         printf("sample3\t-0.31\n");
         printf("...\n");
         printf("Unkown1\t-9999\n");
         printf("Unkown2\t-9999\n");
         printf("Note: make sure that the order of the samples matches those in the .value file\n\n");
         exit(0); 
      } 
   }
   fclose(fp);
   *sample=alloc_sample(cn);
   for (i=0;i<cn;i++) {
     (*sample)[i].name=p_first->name;
     (*sample)[i].value=p_first->value;
     p_sample=p_first->next;
     free(p_first);
     p_first=p_sample;
   }
   printf("%s\t%d\n",fileName,cn);

   return (cn);
}

// Temporary linked-list for data input.
typedef struct sample_info_list_ {
   float *value;
   char *symbol;
   int count;
   struct sample_info_list_ *next;
} sample_info_list;

// Read data into a temporary linked-list, then transfer into a flat array.
int read_value(Data **data,int numSample,const char *fileName,Samples *sample) {

   FILE *fp;
   int numVariable,len,i,cn;
   // buffer should be  enough for 'numSample' values per line.
   char buffer[12*numSample+256],*tok;
   sample_info_list *p_first=NULL, *p_last, *p_data;

   fp=fopen(fileName,"r");
   if (!fp) { perror(fileName); exit(1); }

   for (cn=0; !feof(fp); cn++) {
      if (!fgets(buffer,sizeof(buffer),fp)) break;
      len=strlen(buffer);
      buffer[len-1]='\0';
      p_data=malloc(sizeof(*p_data));
      if (!p_data) { printf("Error allocating data!\n"); exit(1); }
      if (!p_first) p_first=p_data;
      else p_last->next=p_data;
      p_last=p_data;
      tok=strtok(buffer,"\t");
      p_data->symbol=alloc_strdup(tok);
      p_data->value=alloc_float(numSample);
      for (i=0; i<numSample; i++) {
         tok=strtok(NULL,"\t");
         p_data->value[i]=atof(tok);
      }
   }
   fclose(fp);

   numVariable=cn;

   *data=alloc_data(numVariable);
   cn=0;
   do {
     if (p_first->symbol) {
       (*data)[cn].symbol=p_first->symbol;
       (*data)[cn].value=p_first->value;
       (*data)[cn].count=0;
       cn++;
     }
     p_data=p_first->next;
     free(p_first);
     p_first=p_data;
     
   } while (p_first);

   printf("%s\t%d\n",fileName,numVariable);
   return (numVariable);
}

