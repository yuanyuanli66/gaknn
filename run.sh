## example of a run script.

populationSize=2000                                     # population size
numGeneration=1000                                      # number of generations
numCycle=100                                            # number of independent cycles/repeats
thread=50                                               # number of threads used
propTest=0.25                                           # proportion of samples set aside as the test samples
chromosome=30                                           # chromosome length (d), d=30
knn=5                                                   # k-nearest neghbors, k=5

outInfo=info.txt                                        # output file name - general
outChr=top_chromosome.txt                               # output file name - top chromosomes, one per cycle/repeat
outPred=cumulative_prediction.txt                       # output file name - prediction 
outCount=selection_count.txt                            # output file name - gene selection and count
outAccuracy=accuracy.txt                                # output file namn - training and testing rmsd, pearson, spearman correlation

dat_folder="./Data"                                     # folder path containing the data
exe_folder="./Code"                                     # folder path for the code
data_list=("1372_trametinib")                           # name of the data files
#data_list=("nci60_354462_Hypothemycin")                # name of the data files
                                                        # in this example, there are two sets of data in the folder 
                                                        # file name: 1372_trametinib with extension .ic50 and .value
                                                        # file namen: test with extension .ic50 and .value
                                                        # the extension can be changed but must match the specified names below
data_list_len=${#data_list[*]}                          # number of datasets to be run (1 in this example)

for p in `seq 0 1 $((${data_list_len}-1))`; do
   file_name=${data_list[$p]}
   mkdir ${data_list[$p]}
   out_folder="./${data_list[$p]}"
   $exe_folder/gaknn \
    -classFile $dat_folder/$file_name.ic50 \
    -dataFile  $dat_folder/$file_name.value \
    -knn $knn \
    -chromosome  $chromosome \
    -propTest    $propTest \
    -thread      $thread \
    -popSize     $populationSize \
    -numGen      $numGeneration \
    -numCycle    $numCycle \
    -outInfo     $out_folder/$outInfo \
    -outChr      $out_folder/$outChr  \
    -outPred     $out_folder/$outPred \
    -outCount    $out_folder/$outCount \
    -outAccuracy $out_folder/$outAccuracy
done

