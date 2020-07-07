# Gene Selection and Sample Prediction using a Genetic Algorithm and K-Nearest Neighbors Algorithm (GA/KNN)


**Leping Li** <br>
**National Institute of Environmental Health Sciences, NIH** <br>
**Durham, North Carolina 27709** <br>
**Email: [li3@niehs.nih.gov](li3@niehs.nih.gov)** <br>


This is an abbreviated documentation of the GA/KNN algorithm. This version of the software (gaknn) can only be used to predict the values of new data points based on how closely they resemble the points in the training set using the k-nearest neighbors (knn) classification method. It cannot perform sample classification—e.g. normal vs tumor—although it can easily be modified for that purpose. For details of the GA/KNN algorithm, see [here](https://academic.oup.com/bioinformatics/article/17/12/1131/225290).

To compile, go to the **Code** folder, type '**make clean**' and then '**make**'. This should generate the executable '**gaknn**'.

gaknn requires two input data files:

* one contains sample names and outcome labels (outcome variable)
* one contains expression data (predictors - matrix)

In the **Data** folder, you can find an example of the input files: *1372\_trametinib.ic50* and *1372\_trametinib.value*. In this example, ***1372_trametinib*** is the file name. Note that the two files have the same file name but different extensions. You can change the extensions to whatever you like but they must match the extensions specified in the run script (example **run.sh** provided) or the file names in the command line (see an example, below). 


* *1372\_trametinib.ic50* contains cell line names and the respective IC<sub>50</sub> values for a drug
* *1372\_trametinib.value* contains the gene expression data for the cell lines

Each row of the expression data (\*.value) corresponds to a gene. Each column corresponds to the expression values of the genes in a sample (e.g., cell line). The number of columns in the ***\.value*** file must be equal to the number of samples in the ***\.ic50*** file + 1 (gene name column). The orders for which the samples appear in the ***\.value*** file must match those in the ***\.ic50*** file.


As mentioned above, the ***\.ic50*** file has two columns—sample name and outcome (e.g., ln(IC<sub>50</sub>)). In both datasets provided, there are no unknown samples whose values need to be predicted. For those datasets, the **gaknn** algorithm will simply divide the samples randomly into a training and testing set and use the training data to select a set of genes (a "chromosome") whose expression data are most predictive of the IC<sub>50</sub> values of samples in the training set. The identified genes are subsequently used to predict the IC<sub>50</sub> values of the test samples. This process is repeated multiple times (e.g., cycle=100 times). If you have independent test set samples whose IC<sub>50</sub> values need to be predicted, you can specify those samples by assigning their IC<sub>50</sub> values as **-9999** or **NA** in the ***\.ic50*** file. Of course, you would need to have the corresponding gene expression data appended to the ***\.value*** file. If the numbers of samples in both the ***\.ic50*** and ***.value*** do not match, the software will not run. 


A few additional points worth of mentioning:

* The algorithm does not perform any data normalization or standardization. If this needs to be done, it must be carried out outside the algorithm. 

* This version of the **gaknn** algorithm only minimizes the objective function, e.g., the squared sum of the differences between the predicted and observed IC<sub>50</sub> values. If you need to _maximize_ the objective function, you can modify the **selection.c** subroutine accordingly.


| Basic parameters  | |
| ------------- | ----- |
| populationSize=5000      | # population size |
| numGeneration=1000 |  # number of generations|
| numCycle=100 | # number of independent cycles|
| thread=60 | # number of threads |
| propTest=0.25 | # proportion of samples with known IC<sub>50</sub> values set aside as the test samples|
| chromosome=30 | # chromosome length (d), d=30 |
| knn=5 | # k-nearest neighbors, k=5 |



| Output file names  | |
| ------------- | ----- |
| outInfo=info.txt | # output file name - general |
| outChr=top_chromosome.txt | # output file name - top chromosomes, one per cycle/repeat |
| outPred=cumulative_prediction.txt | # output file name - prediction  |
| outCount=selection_count.txt | # output file name - gene selection and count |
| outAccuracy=accuracy.txt | # output file name - training and testing rmsd, Pearson, Spearman correlation |

You can simply run the algorithm using default settings as follows (assuming the executable **gaknn** is in ./**Code** and the data files are in ./**Data**) folders.

./Code/gaknn -classFile ./Data/1372\_trametinib.ic50 -dataFile ./Data/1372_trametinib.value

For those who are familiar with shell script, you can modify the included '**run.sh**' file accordingly.

For questions and comments, send an email to Leping Li at [li3@niehs.nih.gov](li3@niehs.nih.gov)








