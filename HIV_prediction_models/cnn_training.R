###Extracted from Steiner et al. impplementation (https://github.com/maggiesteiner/HIV_DeepLearning.git)
### CNN Run Script ###

#required packages
require(keras3)
require(caret)
require(iml)
require(plyr)
set.seed(1234)

require(reticulate)
use_condaenv("use_tensorflow")

#functions

get_model<-function(bp,nunits, vocab_size = 31){
  #create model
  model<-keras_model_sequential() %>%
    layer_embedding(input_dim=bp,output_dim=128)%>%#,input_length = bp) %>%
    layer_conv_1d(filters=32,kernel_size=9,activation="relu") %>%
    layer_max_pooling_1d(pool_size = 5)%>%
    layer_conv_1d(filters=32,kernel_size=9,activation="relu") %>%
    layer_flatten() %>%
    layer_dense(units=1,activation="sigmoid")

  #and we compile the model ##modified for keras3
  model$compile(
    optimizer = optimizer_rmsprop(),
    loss = "binary_crossentropy",
    metrics = list("accuracy")
  )
  summary(model)
  return(model)
}

#list of input files
input_list = c("atv.fasta","ddi.fasta","etr.fasta","rpv.fasta","drv.fasta","tpv.fasta",
               "lpv.fasta","sqv.fasta","nfv.fasta","fpv.fasta","idv.fasta","d4t.fasta",
               "abc.fasta","azt.fasta","efv.fasta","nvp.fasta","tdf.fasta","3tc.fasta")

drug_list = c("ATV","DDI","ETR","RPV","DRV","TPV","LPV","SQV","NFV","FPV","IDV","D4T","ABC","AZT","EFV","NVP","TDF","3TC") #changed by MSM

#iterate through input files
for (i in 1:length(input_list)){
  drug = drug_list[i]

  infile_train=paste0("../datasets/fasta",input_list[i])
  infile<-gsub(".fasta","", infile_train)
  if (i %in% c(1,5,6,7,8,9,10,11)) genecode = "PR"
  else genecode="RT"

  if (i %in% c(1,5,6,7,8,9,10,11)) dataset = "PI"
  else if (i %in% c(3,4,15,16)) dataset = "NNRTI"
  else dataset = "NRTI"
  
  outfile=paste(infile,".cnn.output.txt",sep="")
  message(genecode)
  #begin writing to outfile
  sink(outfile)
  
  ### data input and formatting ###
  
  if(genecode=="RT"){
    bp<-as.integer(240)
    nunits<-99
  }
  if(genecode=="PR"){
    bp<-as.integer(99)
    nunits<-33
  }
  fasta_train<-read.table(infile_train,stringsAsFactors = FALSE,comment.char = "")
  data_train<-array(dim=c(nrow(fasta_train),2))
  seqid_train<-array(dim=c(nrow(fasta_train)/2,1))

  for(i in 1:nrow(fasta_train)){
    if(i%%2==0){ #even (sequences)
      data_train[i-1,2]<-fasta_train[i,]
    }
    else #odd
    {
      data_train[i,1]<-fasta_train[i,]
    }
  }#end for
  evens_train<-seq(0,nrow(fasta_train),by=2)
  data_train<-data_train[-evens_train,]
  for(i in 1:nrow(data_train)){ #strip IDs
    seqid_train[i]<-substr(data_train[i,1],2,nchar(data_train[i,1])-2)
    seqid_train[i]<-substr(seqid_train[i],1,nchar(seqid_train[i]))
    data_train[i,1]<-substr(data_train[i,1],nchar(data_train[i,1]),nchar(data_train[i,1]))
  }#end for


  data_labels_train<-data_train[,1]
  data_seqs_train<-data_train[,2]
  
  data_labels_train<-as.numeric(data_labels_train)
  
  #convert characters to integers
  seqs_num_train<-array(dim=c(length(data_seqs_train),bp))
  for(i in 1:length(data_seqs_train)){
    z<-data_seqs_train[i]
    seq<-unlist(strsplit(z,""))
    for(k in 1:length(seq)){
      seq[k]<-switch(seq[k],
                     "A"=1,
                     "a"=1,
                     "B"=2,
                     "b"=2,
                     "C"=3,
                     "c"=3,
                     "D"=4,
                     "d"=4,
                     "E"=5,
                     "e"=5,
                     "F"=6,
                     "f"=6,
                     "G"=7,
                     "g"=7,
                     "H"=8,
                     "h"=8,
                     "I"=9,
                     "i"=9,
                     "J"=10,
                     "j"=10,
                     "K"=11,
                     "k"=11,
                     "L"=12,
                     "l"=12,
                     "M"=13,
                     "m"=13,
                     "N"=14,
                     "n"=14,
                     "O"=15,
                     "o"=15,
                     "P"=16,
                     "p"=16,
                     "Q"=17,
                     "q"=17,
                     "R"=18,
                     "r"=18,
                     "S"=19,
                     "s"=19,
                     "T"=20,
                     "t"=20,
                     "U"=21,
                     "u"=21,
                     "V"=22,
                     "v"=22,
                     "W"=23,
                     "w"=23,
                     "X"=24,
                     "x"=24,
                     "Y"=25,
                     "y"=25,
                     "Z"=26,
                     "z"=26,
                     "."=27,
                     "#"=28,
                     "~"=29,
                     "*"=30,
                     0
      )
    }
    seqs_num_train[i,]<-as.integer(seq)
  }
  
  #convert data into list
  data_list_train<-c()
  for(i in 1:nrow(seqs_num_train)){
    seqi<-seqs_num_train[i,]
    data_list_train[[i]]<-seqi
  }

  #pad sequences
  data_f_train <- pad_sequences(
    sequences = data_list_train,
    padding = "post",
    maxlen = bp
  )

  #set number of folds
  k<-5 ##we use the 5-fold cross validation

  #k-fold cross validation
  starttime=Sys.time()
  for(j in seq(0,k-1,1)){

    ##We separate the data into training and testing sets
    fold_data <- read.table(file = paste0("../datasets/", dataset,"/", dataset, "_", drug, "_5folds.tsv"), header = TRUE, sep = "\t" ,stringsAsFactors = FALSE)
    train_fold_seqs = fold_data[fold_data$fold != j,]$SeqID

    training_data = data_f_train[which(seqid_train %in% train_fold_seqs),]
    training_labels = data_labels_train[which(seqid_train %in% train_fold_seqs)]

    #class weights
    zero=length(which(training_labels==0))
    one=length(which(training_labels==1))
    weight_0 = 1 
    weight_1 = zero/one 

    model<-get_model(bp,nunits)

    #train model
    model %>% fit(
      x=training_data,
      y=training_labels,
      epochs=500,
      batch_size=64,
      class_weight = list(`0` = weight_0, `1` = weight_1)
    )
    summary(model)
    
    # We export the model
    model$save(paste("cnn_models/", infile,"_model_fold_", j, ".keras",sep=""))
    message(paste(infile, "model saved!", sep=""))

  }#end of iteration through folds
  
  print(validation_scores)
  endtime = Sys.time()
  print("Runtime:")
  print(endtime-starttime)

}
