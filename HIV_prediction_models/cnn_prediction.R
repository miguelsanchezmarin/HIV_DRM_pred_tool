###Extracted from Steiner et al. impplementation (https://github.com/maggiesteiner/HIV_DeepLearning.git)
###We import the CNN models and predict the test sets
### CNN Run Script ###

#we set the python 3.10 environment
require(reticulate)
use_condaenv("use_tensorflow")

#required packages
require(keras3)

drugs = c('EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV')
NNRTI = c('EFV', 'NVP', 'ETR', 'RPV')
NRTI = c('3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF')
PI = c('FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV')

for(d in drugs){
    for(f in seq(0,4,1)){#we iterate over the 5 folds
        #we load the model
        model <- keras3::load_model(paste0("cnn_models/",tolower(d), "_model_fold_", f, ".keras"))

        #if the drug is in the NNRTI or NRTI class, we set the bp to 240
        if(d %in% NNRTI){
            bp<-as.integer(240)
            nunits<-80
        }

        if(d %in% NRTI){
            bp<-as.integer(240)
            nunits<-80
        }

        if(d %in% PI){
            bp<-as.integer(99)
            nunits<-33
        }
        

        #we load the FASTA test set
        infile_test <- paste0("../datasets/fasta/",d, "_test_fold_", f, ".fasta")
        fasta_test <- read.table(infile_test,stringsAsFactors = FALSE,comment.char = "")
        data_test<-array(dim=c(nrow(fasta_test),2))
        seqid_test<-array(dim=c(nrow(fasta_test),1))
        for(i in 1:nrow(fasta_test)){
        if(i%%2==0){ #even (sequences)
            data_test[i-1,2]<-fasta_test[i,]
        }
        else #odd
        {
        data_test[i,1]<-fasta_test[i,]
        }
        }#end for
        evens_test<-seq(0,nrow(fasta_test),by=2)
        data_test<-data_test[-evens_test,]
        for(i in 1:nrow(data_test)){ #strip IDs
            #we store the SeqID (removing the first and the two last elements)
            seqid_test[i]<-substr(data_test[i,1],2,nchar(data_test[i,1])-2)
            seqid_test[i]<-substr(seqid_test[i],1,nchar(seqid_test[i]))
            data_test[i,1]<-substr(data_test[i,1],nchar(data_test[i,1]),nchar(data_test[i,1]))
        }#end for
  
        data_labels_test<-data_test[,1]
        data_seqs_test<-data_test[,2]
        data_labels_test<-as.numeric(data_labels_test)
    
        seqs_num_test<-array(dim=c(length(data_seqs_test),bp))
        for(i in 1:length(data_seqs_test)){
            z<-data_seqs_test[i]
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

            seqs_num_test[i,]<-as.integer(seq)
        }
        data_list_test<-c()
        for(i in 1:nrow(seqs_num_test)){
            seqi<-seqs_num_test[i,]
            data_list_test[[i]]<-seqi
        }
        data_f_test <- pad_sequences(
        sequences = data_list_test,
        padding = "post",
        maxlen = bp
        )

        #we predict the test set
        preds_prob <- model %>% predict(data_f_test)
        preds_test <- ifelse(preds_prob > 0.5, 1, 0)

        #we save a table with the SeqID and the prediction
        preds_test_df <- data.frame(SeqID = seqid_test, CNN_Prediction = preds_test, CNN_Prob = preds_prob, Real_Prediction = data_labels_test)
        write.table(preds_test_df, file = paste0("../method_predictions/cnn/cnn_", d, "_fold_", f, "predictions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
}