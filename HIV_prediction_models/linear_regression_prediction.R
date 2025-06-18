###ADAPTED FROM https://hivdb.stanford.edu/pages/genopheno.dataset.html. Additions for the project are marked as 'added by MSM'.
# Original script from:
# author: Haley Hedlin
# email: hedlin@stanford.edu
# last updated: September 24, 2014


# To run this function:
# 1. Open R (Install the caret R package if not already installed)
# 2. Open the DRMcv.R file and run the code contained within.
# 3. Type DRMcv() at the R prompt
# 4. Two files will be output into your R working directory

### arguments to DRMcv

# dataset specifies which dataset the function is to be run on
# must be either "PR", "NRTI", or "NNRTI"

# muts.in is a character vector of the mutations to be included as
# independent variables in the OLS
# each entry in the vector is 3 characters long
  # a letter representing the mutation in the first character
  # followed by two numbers representing the position

# drug is a string indicating the drug to be used for the dependent variable
# must be equal to one of the following: "FPV", "ATV", "IDV", "LPV", "NFV", 
# "RTV", "SQV", "TPV", "DRV"

# min.muts is the minimum number of sequences that a mutation must appear in.
# If a mutation appears in too few sequences, it is removed from the model.

# nfold is the number of folds in the cross-validatioin (CV)
# nrep is the number of times to repeat the CV

# confusion controls whether the confusion matrix should be output

# lars controls whether LARS estimates are output


### output from DRMcv.cleanin

# DRMcv outputs up to four files 
# the files will write to R's working directory
# type getwd() in R to see the working directory path

# the first file contains a matrix containing the OLS estimates
# the estimated coefficients for the input mutations are in the first column
# and the SEs of the estimates are in the second column

# the second file contains the mean square errors (MSEs) estimated from the CV folds

# if confusion=TRUE, the third file will contain a confusion matrix 
# with the actual class in the rows and the predicted class in the columns

# if lars=TRUE, the last file will contain the coefficients from LARS 


DRMcv <- function(dataset="PI", drug="LPV", min.muts=10, nfold=5, nrep=10,
                  muts.in=c("47A", "84A", "50V", "76V", "82A", "82F", "84V",
                            "84C", "82S", "82T", "82M", "32I", "47V", "54M", "54L",
                            "54V", "90M", "54A", "54S", "54T", "46I", "46L", "48V",
                            "48M", "24I", "82C", "33F", "10F", "73S", "73T", "73C",
                            "73A", "11I", "11L", "89V", "20T", "53L", "88S", "50L",
                            "24F", "30N", "43T", "46V", "58E", "83D", "88T", "85V"),
                       add_combinations=FALSE, #added by MSM 
                       confusion=FALSE,lars=FALSE){
  
  require(caret)
  require(dplyr)
  
  # check that arguments are entered correctly
  if(!is.character(muts.in)){
    stop('The argument "muts.in" must be a character vector.')
  }
  if(!is.character(drug) | length(drug)!=1){
    stop('The argument "drug" must be a character vector of length 1.')
  }
  
  if(!is.character(dataset) | !dataset%in%c("PI","NRTI","NNRTI", "INI")){
    stop('The argument "dataset" must be a character vector equal to "PI", 
         "NRTI", "NNRTI", or "INI".') ##changed by MSM
  }
  
  drg <- drug 
  if(drug=="3TC") drug <- "X3TC"
  
  muts.in.conv <- convert.muts(muts.in)
  check.muts(muts.in.conv)
  
  # get the amino acids and positions for the mutations to be included in the model
  mut <- ifelse(nchar(muts.in)==3,toupper(substr(muts.in,3,3)),
                toupper(substr(muts.in,4,4)))
  ps <- suppressWarnings(ifelse(nchar(muts.in)==3,as.numeric(substr(muts.in,1,2)),
                                as.numeric(substr(muts.in,1,3)))) 



  for(k in seq(0,4,1)){ #we loop over the 5 folds
    data <- read.table(file = paste0("../datasets/", dataset, "_", drg, "_5folds.tsv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    ##we use as test set the rows with fold==k
    test_seq <- data[data$fold==k,]
    ##we use as train set the rows with fold!=k
    train_seq <- data[data$fold!=k,]

    ## automatically read in the data using url
    if(dataset=="PI"){ 
      dat <- read.table("../datasets/PI_dataset.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

      #we keep in the dat for the training the sequences with SeqID in train_seq SeqID in a dataframe format
      train <- dat[dat$SeqID %in% train_seq$SeqID,]
      test <- dat[dat$SeqID %in% test_seq$SeqID,]
      
      posu <- train[,10:108]
      posu_test <- test[,10:108]
    }
    if(dataset=="NRTI"){
      dat <- read.table("../datasets/NRTI_dataset.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

      #we keep in the dat for the training the sequences with SeqID in train_seq SeqID
      train <- dat[dat$SeqID %in% train_seq$SeqID,]
      test <- dat[dat$SeqID %in% test_seq$SeqID,]
 
      posu <- train[,8:247]
      posu_test <- test[,8:247]
    }
    if(dataset=="NNRTI"){ 
      dat <- dat <- read.table("../datasets/NNRTI_dataset.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
      #we keep in the dat for the training the sequences with SeqID in train_seq SeqID
      train <- dat[dat$SeqID %in% train_seq$SeqID,]
      test <- dat[dat$SeqID %in% test_seq$SeqID,]
 
      posu <- train[,6:245]
      posu_test <- test[,6:245]
    }
    if(dataset=="INI"){ #added by MSM
      dat <- read.table("../datasets/INI_dataset.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

      #we keep in the dat for the training the sequences with SeqID in train_seq SeqID
      train <- dat[dat$SeqID %in% train_seq$SeqID,]
      test <- dat[dat$SeqID %in% test_seq$SeqID,]
 
      posu <- train[,7:294]
      posu_test <- test[,7:294]
    }

    # construct design matrix for OLS
    X <- buildX(posu, mut, ps)
    X_test <- buildX(posu_test, mut, ps)
  
    if(add_combinations==TRUE){ ##we add the combinations described in StanfordHIVDB as variables #added by MSM
      X <- add_X_combinations(X, train, dataset)
      X_test <- add_X_combinations(X_test, test, dataset)
    }

    ##The same but for the train
    drugcol <- which(names(train)==drug)
    Y <- as.numeric(train[,drugcol])  # absolute measure
    Ylog10 <- log10(Y)
    Ylog10[is.infinite(Ylog10)] <- NA #changed by MSM
    df.log <- data.frame(Y=Ylog10, X=X)

    #and for the test
    drugcol_test <- which(names(test)==drug)
    Y_test <- as.numeric(test[,drugcol_test])  # absolute measure
    Ylog10_test <- log10(Y_test)
    Ylog10_test[is.infinite(Ylog10_test)] <- NA #changed by MSM
    df.log_test <- data.frame(Y=Ylog10_test, X=X_test)

    #we get the row indexes of the na in the test
    seqID_test <- test[,"SeqID"]
    df.log.cc <- na.omit(df.log) ##modified as we already looked for this beforehand
    df.log_test.cc <- df.log_test

    # remove mutations that are rare
    rare.muts <- which(colSums(df.log.cc[,-1])<min.muts)
    col_keys <- colnames(df.log.cc[,-1])

    if(length(rare.muts)>0){
      message(paste0(substring(col_keys[rare.muts],3),
                     " excluded from the model because it appears in fewer than ",
                     min.muts," sequences.\n"))
      df.log.cc <- df.log.cc[,-(rare.muts+1)]
      df.log_test.cc <- df.log_test.cc[,-(rare.muts+1)]
    }

    ### fit the models 
    fit <- lm(Y~., data=df.log.cc, na.action=na.exclude) #changed by MSM
  
    CVout <- train(df.log.cc[,-1], df.log.cc$Y, method='lm', 
                   trControl=trainControl(method='repeatedcv',number=nfold,repeats=nrep)) ##modified by MSM
  
    if(confusion==TRUE){
      ### code for the confusion matrix
      cutoffmat <- matrix(NA, nrow=22, ncol=2)
      rownames(cutoffmat) <- c("FPV","ATV","IDV","LPV","NFV","SQV","TPV","DRV",
                               "X3TC","ABC","AZT","D4T","DDI","TDF",
                               "EFV","NVP","ETR","RPV",
                               "EVG", "RAL", "DTG", "BIC") #added by MSM
      colnames(cutoffmat) <- c("lower","upper")
      cutoffmat[1,] <- c(3,15) # FPV 
      cutoffmat[2,] <- c(3,15) # ATV
      cutoffmat[3,] <- c(3,15) # IDV
      cutoffmat[4,] <- c(9,55) # LPV
      cutoffmat[5,] <- c(3,6) # NFV
      cutoffmat[6,] <- c(3,15) # SQV
      cutoffmat[7,] <- c(2,8) # TPV
      cutoffmat[8,] <- c(10,90) # DRV
      cutoffmat[9,] <- c(3,25) # X3TC
      cutoffmat[10,] <- c(2,6) # ABC 
      cutoffmat[11,] <- c(3,15) # AZT
      cutoffmat[12,] <- c(1.5,3) # D4T
      cutoffmat[13,] <- c(1.5,3) # DDI
      cutoffmat[14,] <- c(1.5,3) # TDF
      cutoffmat[15,] <- c(3,10) # EFV
      cutoffmat[16,] <- c(3,10) # NVP
      cutoffmat[17,] <- c(3,10) # ETR
      cutoffmat[18,] <- c(3,10) # RPV
      #added by MSM
      cutoffmat[19,] <- c(4,15) # EVG
      cutoffmat[20,] <- c(4,15) # RAL
      cutoffmat[21,] <- c(3.5,15) # DTG
      cutoffmat[22,] <- c(3.5,15) # BIC

    
      cutoff <- cutoffmat[which(rownames(cutoffmat)==drug),]
  
      # predicted and actual categories ##modified for prediction in test
      predicted <- cut(10^predict(fit, newdata = df.log_test.cc),c(0,cutoff,Inf),labels=FALSE)
      actual <- cut(10^df.log_test.cc$Y,c(0,cutoff,Inf),labels=FALSE)

      #we store the SeqID, predicted value and label in a dataframe
      df_confusion <- data.frame(SeqID = seqID_test,RF=10^predict(fit, newdata = df.log_test.cc))

      #we create with this a confusion matrix for labels 1, 2 and 3
      #we add null values for the missing labels
      predicted <- factor(predicted, levels=1:3)
      actual <- factor(actual, levels=1:3)

      conftab <- table(predicted,actual)
      rownames(conftab) <- colnames(conftab) <- c("susceptible",
                                                  "intermediate-level resistant",
                                                  "high-level resistant")
    
    }
  
  
    if(lars==TRUE){
      require(glmnet)
    
      larsfit <- cv.glmnet(as.matrix(df.log.cc[,-1]), df.log.cc$Y, nfolds=5, na.action=na.exclude)##changed by MSM
      larscoef <- coef(larsfit,s="lambda.min")
    
      #added by MSM to predict also with LARS
      if(confusion==TRUE){
        predicted_lars <- cut(10^predict(larsfit, newx=as.matrix(df.log_test.cc[,-1]), s="lambda.min"),c(0,cutoff,Inf),labels=FALSE)
        predicted_lars <- factor(predicted_lars, levels=1:3)
        conftab_lars <- table(predicted_lars,actual)
        rownames(conftab_lars) <- colnames(conftab_lars) <- c("susceptible",
                                                              "intermediate-level resistant",
                                                              "high-level resistant")

        df_confusion_lars <- data.frame(SeqID = seqID_test,RF=10^predict(larsfit, newx=as.matrix(df.log_test.cc[,-1]), s="lambda.min"))
      }
    }
  
    # output model coefficients and SEs, and the MSE  
    if(add_combinations==TRUE) drug <- paste0("_I_",drug) #added by MSM
    # write.table(summary(fit)$coefficients[,1:2],file=paste0("OLS_",drug,"_tsm_fold_", k, ".txt"))
    # write.table(CVout$resample$RMSE^2,file=paste0("CVmse_",drug,"_tsm_fold_", k, ".txt"))
    # write.table(CVout$resample$Rsquared,file=paste0("CVrsq_",drug,"_tsm_fold_", k, ".txt"))
    # if(confusion==TRUE) write.table(conftab,file=paste0("Confusion_",drug,"_tsm_fold_", k, ".txt"))
    # if(lars==TRUE) write.table(cbind(rownames(larscoef),matrix(larscoef)),file=paste0("LARS_",drug,"_tsm_fold_", k, ".txt"))
    # if(confusion==TRUE & lars==TRUE) write.table(conftab_lars,file=paste0("Confusion_LARS_",drug,"_tsm_fold_", k, ".txt")) ##added by MSM
    # write.table(conftab,file=paste0("Confusion_",drug,"_tsm_fold_", k, ".txt"))
    write.table(df_confusion, file=paste0("../method_predictions/linear_regression/LSR",drug,"_fold_", k, "_predictions.txt")) ##added by MSM
    if(lars==TRUE) write.table(df_confusion_lars, file=paste0("../method_predictions/linear_regression/LARS",drug,"_fold_", k, "_predictions.txt")) ##added by MSM

    #we reset the drug string
    if(add_combinations==TRUE) drug <- substr(drug, 1, nchar(drug)-13) #added by MSM
  }

}
#### helper functions

# function to check that the mutations have been entered correctly
# if the letters are entered as lower case, they are converted to upper case

# convert insertions or deletions to # or ~
convert.muts <- function(muts.in){
  muts.in5 <- which(nchar(muts.in) == 5)
  muts.in6 <- which(nchar(muts.in) == 6)
  for(mi in muts.in5){
    postmp <- substr(muts.in[mi],1,2)
    if(substr(muts.in[mi],3,5) == "ins") muts.in[mi] <- paste0(postmp,"#")
    if(substr(muts.in[mi],3,5) == "del") muts.in[mi] <- paste0(postmp,"~")   
  }
  for(mi in muts.in6){
    postmp <- substr(muts.in[mi],1,3)
    if(substr(muts.in[mi],4,6) == "ins") muts.in[mi] <- paste0(postmp,"#")
    if(substr(muts.in[mi],4,6) == "del") muts.in[mi] <- paste0(postmp,"~")   
  }
  return(muts.in)
}

check.muts <- function(muts.in){
  # all entries should be nchar=3
  if(any(nchar(muts.in) < 3) | any(nchar(muts.in) > 4))
    stop('All entries in argument "muts.in" should be between 3 and 4 characters long.')
  
  muts.in3 <- muts.in[nchar(muts.in) == 3]
  muts.in4 <- muts.in[nchar(muts.in) == 4]
  
  # all should have numbers first two or three characters and a letter for the last 
  if(!all(toupper(substr(muts.in3,3,3))%in%c(LETTERS,"#","~")))
    stop('All entries in argument "muts.in" must have a letter, #, or ~ in the last 
         character.')
  if(any(is.na(as.numeric(substr(muts.in3,1,2)))))
    stop('All entries in argument "muts.in" must begin in two or three digits.')
  
  if(!all(toupper(substr(muts.in4,4,4))%in%c(LETTERS,"#","~")))
    stop('All entries in argument "muts.in" must have a letter, #, or ~ in the last 
         character.')
  if(any(is.na(as.numeric(substr(muts.in4,1,3)))))
    stop('All entries in argument "muts.in" must begin in two or three digits.')
}



# function to create the design matrix X with the input mutations/positions
buildX <- function(dat, mut, ps){
  X <- matrix(NA, nrow=nrow(dat), ncol=length(mut))

  # loop through all positions
  for(p in unique(ps)){
    p1 <- substr(dat[,p],1,1)  # first mutation at this position
    p2 <- substr(dat[,p],2,2)
    for(ind in which(ps==p)){
      X[,ind] <-  as.numeric(p1==as.character(mut[ind]) | 
                               p2==as.character(mut[ind]))  
    }
  }  
  colnames(X) <- paste0(ps,mut)  
  return(X)
}

##function to add combinations of mutations as columns to the mutation matrix (added by MSM)
add_X_combinations <- function(X, dat, dataset="PI"){
  print(colnames(X))
  colnames_list <- paste(colnames(X), collapse = ' ')
  colnames_list <- paste0(" ", colnames_list, " ")
  library("stringr")
  X_comb <- X
  if(dataset=="PI"){
    ##we add the combinations of mutations as variables (when they appear more than 10 times) for the PI dataset
    ###we extracted the described combinations from StanfordHIVDB on February 12th (https://hivdb.stanford.edu/dr-summary/mut-scores/PI/)
    p11 <- substr(dat[,which(colnames(dat)=="P11")],1,1)
    p32 <- substr(dat[,which(colnames(dat)=="P32")],1,1)
    p46 <- substr(dat[,which(colnames(dat)=="P46")],1,1)
    p47 <- substr(dat[,which(colnames(dat)=="P47")],1,1)
    p53 <- substr(dat[,which(colnames(dat)=="P53")],1,1)
    p54 <- substr(dat[,which(colnames(dat)=="P54")],1,1)
    p73 <- substr(dat[,which(colnames(dat)=="P73")],1,1)
    p76 <- substr(dat[,which(colnames(dat)=="P76")],1,1)
    p82 <- substr(dat[,which(colnames(dat)=="P82")],1,1)
    p84 <- substr(dat[,which(colnames(dat)=="P84")],1,1)
    p89 <- substr(dat[,which(colnames(dat)=="P89")],1,1)
    p90 <- substr(dat[,which(colnames(dat)=="P90")],1,1)
    #we add a column (11IL+32I) that is 1 if 11I and 32I are present or 11L and 32I are present
    new_col <- as.numeric(str_detect(p11, "^[IL]$") & str_detect(p32, "^I$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "11IL+32I"
    
    ##Now we repeat this for the rest of the combinations
    ##11IL+54LM
    new_col <- as.numeric(str_detect(p11, "^[IL]$") & str_detect(p54, "^[LM]$")) 
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "11IL+54LM"
    
    ##32I+47AV
    new_col <- as.numeric(str_detect(p32, "^I$") & str_detect(p47, "^[AV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "32I+47AV"

    ##32I+54LM
    new_col <- as.numeric(str_detect(p32, "^I$") & str_detect(p54, "^[LM]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "32I+54LM"

    ##32I+76V
    new_col <- as.numeric(str_detect(p32, "^I$") & str_detect(p76, "^V$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "32I+76V"
    
    ##32I+84V
    new_col <- as.numeric(str_detect(p32, "^I$") & str_detect(p84, "^V$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "32I+84V"
    
    ##32I+89V
    new_col <- as.numeric(str_detect(p32, "^I$") & str_detect(p89, "^V$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "32I+89V"

    ##46IL+84V+90M
    new_col <- as.numeric(str_detect(p46, "^[IL]$") & str_detect(p84, "^V$") & str_detect(p90, "^M$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "46IL+84V+90M"
    
    ##46ILV+82ACFLMST
    new_col <- as.numeric(str_detect(p46, "^[ILV]$") & str_detect(p82, "^[ACFLMST]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "46ILV+82ACFLMST"
    
    ##46ILV+90M
    new_col <- as.numeric(str_detect(p46, "^[ILV]$") & str_detect(p90, "^M$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "46ILV+90M"
    
    ##46ILV+76V
    new_col <- as.numeric(str_detect(p46, "^[ILV]$") & str_detect(p76, "^V$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "46ILV+76V"
    
    ##47AV+54LM
    new_col <- as.numeric(str_detect(p47, "^[AV]$") & str_detect(p54, "^[LM]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "47AV+54LM"
    
    ##47AV+84V
    new_col <- as.numeric(str_detect(p47, "^[AV]$") & str_detect(p84, "^V$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "47AV+84V"
    
    ##53L+90M
    new_col <- as.numeric(str_detect(p53, "^L$") & str_detect(p90, "^M$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "53L+90M"

    ##54ALMSTV+82ACFLMST
    new_col <- as.numeric(str_detect(p54, "^[ALMSTV]$") & str_detect(p82, "^[ACFLMST]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "54ALMSTV+82ACFLMST"
    
    ##54ALMSTV+90M
    new_col <- as.numeric(str_detect(p54, "^[ALMSTV]$") & str_detect(p90, "^M$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "54ALMSTV+90M"
    
    ##54LM+84V
    new_col <- as.numeric(str_detect(p54, "^[LM]$") & str_detect(p84, "^V$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "54LM+84V"
    
    ##54LM+89V
    new_col <- as.numeric(str_detect(p54, "^[LM]$") & str_detect(p89, "^V$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "54LM+89V"
    
    ##73ACSTV+90M
    new_col <- as.numeric(str_detect(p73, "^[ACSTV]$") & str_detect(p90, "^M$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "73ACSTV+90M"
    
    ##82ACFLMST+90M
    new_col <- as.numeric(str_detect(p82, "^[ACFLMST]$") & str_detect(p90, "^M$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "82ACFLMST+90M"
    
  }
  if(dataset=="NRTI"){
    p40 <- substr(dat[,which(colnames(dat)=="P40")],1,1)
    p41 <- substr(dat[,which(colnames(dat)=="P41")],1,1)
    p44 <- substr(dat[,which(colnames(dat)=="P44")],1,1)
    p62 <- substr(dat[,which(colnames(dat)=="P62")],1,1)
    p65 <- substr(dat[,which(colnames(dat)=="P65")],1,1)
    p67 <- substr(dat[,which(colnames(dat)=="P67")],1,1)
    p68 <- substr(dat[,which(colnames(dat)=="P68")],1,1)
    p69 <- substr(dat[,which(colnames(dat)=="P69")],1,1)
    p70 <- substr(dat[,which(colnames(dat)=="P70")],1,1)
    p74 <- substr(dat[,which(colnames(dat)=="P74")],1,1)
    p77 <- substr(dat[,which(colnames(dat)=="P77")],1,1)
    p115 <- substr(dat[,which(colnames(dat)=="P115")],1,1)
    p116 <- substr(dat[,which(colnames(dat)=="P116")],1,1)
    p151 <- substr(dat[,which(colnames(dat)=="P151")],1,1)
    p184 <- substr(dat[,which(colnames(dat)=="P184")],1,1)
    p210 <- substr(dat[,which(colnames(dat)=="P210")],1,1)
    p215 <- substr(dat[,which(colnames(dat)=="P215")],1,1)
    p219 <- substr(dat[,which(colnames(dat)=="P219")],1,1)
    ##40F+41L+210W+215FY
    new_col <- as.numeric(str_detect(p40, "^F$") & str_detect(p41, "^L$") & str_detect(p210, "^W$") & str_detect(p215, "^[FY]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "40F+41L+210W+215FY"
  
    ##41L+215FY
    new_col <- as.numeric(str_detect(p41, "^L$") & str_detect(p215, "^[FY]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "41L+215FY"
    
    ##41L+210W
    new_col <- as.numeric(str_detect(p41, "^L$") & str_detect(p210, "^W$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "41L+210W"
    
    ##41L+210W+215FY
    new_col <- as.numeric(str_detect(p41, "^L$") & str_detect(p210, "^W$") & str_detect(p215, "^[FY]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "41L+210W+215FY"
    
    ##41L+44AD+210W+215FY
    new_col <- as.numeric(str_detect(p41, "^L$") & str_detect(p44, "^[AD]$") & str_detect(p210, "^W$") & str_detect(p215, "^[FY]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "41L+44AD+210W+215FY"
    
    ##41L+67EGNHST+215FY
    new_col <- as.numeric(str_detect(p41, "^L$") & str_detect(p67, "^[EGNHST]$") & str_detect(p215, "^[FY]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "41L+67EGNHST+215FY"
    
    ##41L+184VI+215FY
    new_col <- as.numeric(str_detect(p41, "^L$") & str_detect(p184, "^[VI]$") & str_detect(p215, "^[FY]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "41L+184VI+215FY"
    
    ##41L+215ACDEILNSV
    new_col <- as.numeric(str_detect(p41, "^L$") & str_detect(p215, "^[ACDEILNSV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "41L+215ACDEILNSV"
    
    ##62V+65R
    new_col <- as.numeric(str_detect(p62, "^V$") & str_detect(p65, "^R$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "62V+65R"
    
    ##65RN+151M
    new_col <- as.numeric(str_detect(p65, "^[RN]$") & str_detect(p151, "^M$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "65RN+151M"
    
    ##65R+68NGDR
    new_col <- as.numeric(str_detect(p65, "^R$") & str_detect(p68, "^[NGDR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "65R+68NGDR"
    
    ##67EGNHST+215FY+219ENQRW
    new_col <- as.numeric(str_detect(p67, "^[EGNHST]$") & str_detect(p215, "^[FY]$") & str_detect(p219, "^[ENQRW]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "67EGNHST+215FY+219ENQRW"
    
    ##67EGNHST+70R+219ENQRW
    new_col <- as.numeric(str_detect(p67, "^[EGNHST]$") & str_detect(p70, "^R$") & str_detect(p219, "^[ENQRW]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "67EGNHST+70R+219ENQRW"
    
    ##67EGNHST+70R+184VI+219ENQRW	
    new_col <- as.numeric(str_detect(p67, "^[EGNHST]$") & str_detect(p70, "^R$") & str_detect(p184, "^[VI]$") & str_detect(p219, "^[ENQRW]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "67EGNHST+70R+184VI+219ENQRW"
    
    ##69#+184VI
    new_col <- as.numeric(str_detect(p69, "^#$") & str_detect(p184, "^[VI]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "69#+184VI"
    
    ##70EGNQST+184IV
    new_col <- as.numeric(str_detect(p70, "^[EGNQST]$") & str_detect(p184, "^[IV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "70EGNQST+184IV"
    
    ##74V+184IV
    new_col <- as.numeric(str_detect(p74, "^V$") & str_detect(p184, "^[IV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "74V+184IV"
    
    ##77L+116Y+151ML
    new_col <- as.numeric(str_detect(p77, "^L$") & str_detect(p116, "^Y$") & str_detect(p151, "^M$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "77L+116Y+151ML"
    
    ##115F+184IV
    new_col <- as.numeric(str_detect(p115, "^F$") & str_detect(p184, "^[IV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "115F+184IV"
    
    ##151M+184IV
    new_col <- as.numeric(str_detect(p151, "^M$") & str_detect(p184, "^[IV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "151M+184IV"
    
    ##210W+215FY
    new_col <- as.numeric(str_detect(p210, "^W$") & str_detect(p215, "^[FY]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "210W+215FY"
    
    ##210W+215ACDEILNSV
    new_col <- as.numeric(str_detect(p210, "^W$") & str_detect(p215, "^[ACDEILNSV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "210W+215ACDEILNSV"
  
  }
  if(dataset=="NNRTI"){
    p98 <- substr(dat[,which(colnames(dat)=="P98")],1,1)
    p100 <- substr(dat[,which(colnames(dat)=="P100")],1,1)
    p101 <- substr(dat[,which(colnames(dat)=="P101")],1,1)
    p103 <- substr(dat[,which(colnames(dat)=="P103")],1,1)
    p106 <- substr(dat[,which(colnames(dat)=="P106")],1,1)
    p108 <- substr(dat[,which(colnames(dat)=="P108")],1,1)
    p138 <- substr(dat[,which(colnames(dat)=="P138")],1,1)
    p179 <- substr(dat[,which(colnames(dat)=="P179")],1,1)
    p181 <- substr(dat[,which(colnames(dat)=="P181")],1,1)
    p184 <- substr(dat[,which(colnames(dat)=="P184")],1,1)
    p188 <- substr(dat[,which(colnames(dat)=="P188")],1,1)
    p190 <- substr(dat[,which(colnames(dat)=="P190")],1,1)
    p221 <- substr(dat[,which(colnames(dat)=="P221")],1,1)
    p225 <- substr(dat[,which(colnames(dat)=="P225")],1,1)
    p227 <- substr(dat[,which(colnames(dat)=="P227")],1,1)
    p234 <- substr(dat[,which(colnames(dat)=="P234")],1,1)
    ##98G+106I
    new_col <- as.numeric(str_detect(p98, "^G$") & str_detect(p106, "^I$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "98G+106I"
    
    ##98G+190ASC
    new_col <- as.numeric(str_detect(p98, "^G$") & str_detect(p190, "^[ASC]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "98G+190ASC"
    
    ##98G+181CIV
    new_col <- as.numeric(str_detect(p98, "^G$") & str_detect(p181, "^[CIV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "98G+181CIV"
    
    ##98G+227CL
    new_col <- as.numeric(str_detect(p98, "^G$") & str_detect(p227, "^[CL]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "98G+227CL"
    
    ##100I+103N
    new_col <- as.numeric(str_detect(p100, "^I$") & str_detect(p103, "^N$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "100I+103N"
    
    ##101E+181C
    new_col <- as.numeric(str_detect(p101, "^E$") & str_detect(p181, "^C$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "101E+181C"
    
    ##101E+190A
    new_col <- as.numeric(str_detect(p101, "^E$") & str_detect(p190, "^A$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "101E+190A"
    
    ##101E+190S
    new_col <- as.numeric(str_detect(p101, "^E$") & str_detect(p190, "^S$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "101E+190S"
    
    ##101E+188L
    new_col <- as.numeric(str_detect(p101, "^E$") & str_detect(p188, "^L$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "101E+188L"
    
    ##101E+184I
    new_col <- as.numeric(str_detect(p101, "^E$") & str_detect(p184, "^I$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "101E+184I"
    
    ##103N+181C
    new_col <- as.numeric(str_detect(p103, "^N$") & str_detect(p181, "^C$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "103N+181C"
    
    ##103N+225H
    new_col <- as.numeric(str_detect(p103, "^N$") & str_detect(p225, "^H$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "103N+225H"
    
    ##103R+179D
    new_col <- as.numeric(str_detect(p103, "^R$") & str_detect(p179, "^D$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "103R+179D"
    
    ##106I+181C
    new_col <- as.numeric(str_detect(p106, "^I$") & str_detect(p181, "^C$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "106I+181C"
    
    ##106I+190S
    new_col <- as.numeric(str_detect(p106, "^I$") & str_detect(p190, "^S$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "106I+190S"
    
    ##106I+234I
    new_col <- as.numeric(str_detect(p106, "^I$") & str_detect(p234, "^I$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "106I+234I"
    
    ##106I+221Y
    new_col <- as.numeric(str_detect(p106, "^I$") & str_detect(p221, "^Y$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "106I+221Y"
    
    ##106A+227CL
    new_col <- as.numeric(str_detect(p106, "^A$") & str_detect(p227, "^[CL]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "106A+227CL"
    
    ##108I+181C
    new_col <- as.numeric(str_detect(p108, "^I$") & str_detect(p181, "^C$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "108I+181C"
    
    ##108I+234I
    new_col <- as.numeric(str_detect(p108, "^I$") & str_detect(p234, "^I$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "108I+234I"
    
    ##138K+184I
    new_col <- as.numeric(str_detect(p138, "^K$") & str_detect(p184, "^I$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "138K+184I"
    
    ##179F+181C
    new_col <- as.numeric(str_detect(p179, "^F$") & str_detect(p181, "^C$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "179F+181C"
    
    ##179T+181C
    new_col <- as.numeric(str_detect(p179, "^T$") & str_detect(p181, "^C$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "179T+181C"
    
    ##181CIV+190ACSTV
    new_col <- as.numeric(str_detect(p181, "^[CIV]$") & str_detect(p190, "^[ACSTV]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "181CIV+190ACSTV"
    
    ##181CIV+221Y
    new_col <- as.numeric(str_detect(p181, "^[CIV]$") & str_detect(p221, "^Y$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "181CIV+221Y"
  }
  if(dataset=="INI"){
    p51 <- substr(dat[,which(colnames(dat)=="P51")],1,1)
    p66 <- substr(dat[,which(colnames(dat)=="P66")],1,1)
    p74 <- substr(dat[,which(colnames(dat)=="P74")],1,1)
    p75 <- substr(dat[,which(colnames(dat)=="P75")],1,1)
    p92 <- substr(dat[,which(colnames(dat)=="P92")],1,1)
    p97 <- substr(dat[,which(colnames(dat)=="P97")],1,1)
    p118 <- substr(dat[,which(colnames(dat)=="P118")],1,1)
    p122 <- substr(dat[,which(colnames(dat)=="P122")],1,1)
    p138 <- substr(dat[,which(colnames(dat)=="P138")],1,1)
    p140 <- substr(dat[,which(colnames(dat)=="P140")],1,1)
    p143 <- substr(dat[,which(colnames(dat)=="P143")],1,1)
    p147 <- substr(dat[,which(colnames(dat)=="P147")],1,1)
    p148 <- substr(dat[,which(colnames(dat)=="P148")],1,1)
    p149 <- substr(dat[,which(colnames(dat)=="P149")],1,1)
    p155 <- substr(dat[,which(colnames(dat)=="P155")],1,1)
    p157 <- substr(dat[,which(colnames(dat)=="P157")],1,1)
    p163 <- substr(dat[,which(colnames(dat)=="P163")],1,1)
    p230 <- substr(dat[,which(colnames(dat)=="P230")],1,1)
    p263 <- substr(dat[,which(colnames(dat)=="P263")],1,1)
    ##51Y+263K
    new_col <- as.numeric(str_detect(p51, "^Y$") & str_detect(p263, "^K$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "51Y+263K"

    ##66I+118R
    new_col <- as.numeric(str_detect(p66, "^I$") & str_detect(p118, "^R$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "66I+118R"
    
    ##74FIM+118R
    new_col <- as.numeric(str_detect(p74, "^[FIM]$") & str_detect(p118, "^R$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "74FIM+118R"
    
    ##74FIM+148HKR
    new_col <- as.numeric(str_detect(p74, "^[FIM]$") & str_detect(p148, "^[HKR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "74FIM+148HKR"
    
    ##74FIM+143ACGHRS
    new_col <- as.numeric(str_detect(p74, "^[FIM]$") & str_detect(p143, "^[ACGHRS]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "74FIM+143ACGHRS"
    
    ##75A+118R
    new_col <- as.numeric(str_detect(p75, "^A$") & str_detect(p118, "^R$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "75A+118R"
    
    ##75A+148HKR
    new_col <- as.numeric(str_detect(p75, "^A$") & str_detect(p148, "^[HKR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "75A+148HKR"
    
    ##92Q+155H
    new_col <- as.numeric(str_detect(p92, "^Q$") & str_detect(p155, "^H$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "92Q+155H"
    
    ##97A+118R
    new_col <- as.numeric(str_detect(p97, "^A$") & str_detect(p118, "^R$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "97A+118R"
    
    ##97A+148HKR	
    new_col <- as.numeric(str_detect(p97, "^A$") & str_detect(p148, "^[HKR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "97A+148HKR"
    
    ##97A+143ACGHRS
    new_col <- as.numeric(str_detect(p97, "^A$") & str_detect(p143, "^[ACGHRS]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "97A+143ACGHRS"
    
    ##118R+138AKT
    new_col <- as.numeric(str_detect(p118, "^R$") & str_detect(p138, "^[AKT]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "118R+138AKT"
    
    ##122N+148HKR
    new_col <- as.numeric(str_detect(p122, "^N$") & str_detect(p148, "^[HKR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "122N+148HKR"
    
    ##138AKT+140ACS
    new_col <- as.numeric(str_detect(p138, "^[AKT]$") & str_detect(p140, "^[ACS]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "138AKT+140ACS"
    
    ##138AKT+148HKR
    new_col <- as.numeric(str_detect(p138, "^[AKT]$") & str_detect(p148, "^[HKR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "138AKT+148HKR"
    
    ##138K+263K
    new_col <- as.numeric(str_detect(p138, "^K$") & str_detect(p263, "^K$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "138K+263K"
    
    ##140ACS+148HKR
    new_col <- as.numeric(str_detect(p140, "^[ACS]$") & str_detect(p148, "^[HKR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "140ACS+148HKR"
    
    ##140ACS+148HKR+149A
    new_col <- as.numeric(str_detect(p140, "^[ACS]$") & str_detect(p148, "^[HKR]$") & str_detect(p149, "^A$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "140ACS+148HKR+149A"
    
    ##143ACGHRS+163R
    new_col <- as.numeric(str_detect(p143, "^[ACGHRS]$") & str_detect(p163, "^R$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "143ACGHRS+163R"
    
    ##143ACGHRS+230R
    new_col <- as.numeric(str_detect(p143, "^[ACGHRS]$") & str_detect(p230, "^R$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "143ACGHRS+230R"
    
    ##147G+148HKR
    new_col <- as.numeric(str_detect(p147, "^G$") & str_detect(p148, "^[HKR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "147G+148HKR"
    
    ##147G+155H
    new_col <- as.numeric(str_detect(p147, "^G$") & str_detect(p155, "^H$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "147G+155H"
    
    ##148HKR+155H
    new_col <- as.numeric(str_detect(p148, "^[HKR]$") & str_detect(p155, "^H$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "148HKR+155H"
    
    ##148HKR+163KR
    new_col <- as.numeric(str_detect(p148, "^[HKR]$") & str_detect(p163, "^[KR]$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "148HKR+163KR"
    
    ##155H+263K
    new_col <- as.numeric(str_detect(p155, "^H$") & str_detect(p263, "^K$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "155H+263K"
    
    ##157Q+263K
    new_col <- as.numeric(str_detect(p157, "^Q$") & str_detect(p263, "^K$"))
    X_comb <- cbind(X_comb, new_col)
    colnames(X_comb)[ncol(X_comb)] <- "157Q+263K"
  }
  return(X_comb)
}

##added by MSM
drugs<- c('RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV')
INIs<- c('RAL', 'EVG', 'DTG', 'BIC')
PIs<- c('FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV')
NNRTIs<- c('EFV', 'NVP', 'ETR', 'RPV')
NRTIs<- c('3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF')

##retrieved with available_mutations.py from the prefiltered dataset
###Complete set of mutations
INI_muts<- c('150A', '207E', '216R', '24T', '196P', '215R', '24N', '283R', '103R', '59A', '163K', '119G', '45S', '90S', '113V', '73V', '69G', '211R', '256E', '161T', '17C', '221E', '70R', '112V', '74I', '261S', '49P', '99H', '72V', '79A', '135V', '163Q', '222T', '143G', '205T', '151A', '147G', '125A', '66K', '198D', '270H', '146K', '14R', '112R', '72A', '187K', '134N', '66I', '205S', '118R', '51Y', '31I', '50V', '279G', '284G', '251L', '153F', '253Y', '125P', '140C', '84M', '258N', '193E', '136R', '283G', '98G', '230N', '181L', '26Y', '55S', '119T', '280S', '193S', '158F', '232E', '125M', '208L', '23S', '23V', '160N', '41G', '122I', '167E', '57G', '111T', '281M', '127R', '146P', '234V', '201I', '151I', '140S', '212T', '39N', '68V', '97A', '92G', '179G', '20K', '260L', '219N', '171Q', '176L', '17T', '218I', '72T', '160R', '148K', '45I', '36V', '30S', '112K', '148R', '95P', '99F', '30A', '124Q', '157K', '84V', '97I', '60L', '182V', '161L', '255G', '126L', '60M', '59E', '162V', '134D', '41E', '255R', '97S', '221H', '200M', '92Q', '279N', '219Q', '240R', '216H', '215Q', '221S', '218Q', '114R', '32I', '124A', '220M', '133T', '28I', '146R', '21S', '154L', '185Y', '91T', '254D', '123G', '163R', '254K', '74M', '24D', '39C', '219T', '270E', '45Q', '115H', '227F', '208M', '148N', '18D', '165I', '230H', '263K', '151L', '111R', '27D', '275V', '112I', '155H', '136Q', '101V', '149A', '13D', '215N', '278A', '119A', '50R', '283D', '57N', '154I', '25E', '112A', '56S', '255N', '174A', '163H', '143H', '41N', '128T', '148H', '27G', '160Q', '188R', '119R', '156N', '270N', '70E', '17N', '261F', '213S', '100Y', '124N', '88I', '95K', '153A', '111Q', '269K', '220L', '286N', '230R', '95R', '22I', '68I', '110T', '221G', '234H', '139Y', '21T', '50L', '80S', '163E', '135L', '60V', '138A', '136N', '136T', '72L', '123C', '111A', '125Q', '262F', '37I', '170A', '70D', '268L', '210I', '231K', '122N', '10A', '232N', '84L', '27S', '254Q', '265V', '163D', '203M', '72N', '217V', '35Q', '153Y', '112S', '210L', '50I', '66A', '10D', '91E', '42R', '206S', '24G', '175T', '212S', '126M', '234I', '212A', '106A', '101I', '119P', '288N', '74V', '24A', '11D', '56Y', '143A', '234F', '63I', '71Q', '143C', '143S', '96D', '222K', '195T', '145S', '54I', '171Y', '141V', '113L', '138K', '210S', '146I', '146L', '157Q', '195C', '216N', '193D', '212Y', '124S', '211Q', '125S', '122V', '34R', '221R', '259L', '63M', '42Q', '262L', '253E', '163T', '218S', '255K', '143R', '200T', '97V', '138D', '121Y', '140A', '171L', '140E', '125V', '212Q', '173R', '193R', '91S', '254S', '229E', '212I', '211T', '220V', '45V', '138T', '50T', '126A')
NNRTI_muts<- c('184I', '138A', '188F', '239I', '123Q', '162D', '162H', '107S', '224N', '64H', '195T', '13R', '22N', '207R', '35E', '121H', '200S', '227L', '109I', '35A', '202V', '211K', '190C', '135P', '215I', '68G', '105L', '82R', '107A', '162Y', '69N', '212M', '36G', '227C', '240E', '166R', '181V', '31R', '219T', '74V', '236R', '136E', '88G', '192N', '225H', '195L', '158T', '41I', '35R', '211T', '164I', '204D', '123E', '200G', '102M', '178V', '228N', '173Q', '211D', '142K', '160L', '103N', '179D', '207H', '178R', '20R', '80I', '69E', '101A', '88C', '138Q', '159V', '136S', '67G', '32N', '204V', '215V', '142L', '100I', '200K', '162C', '108I', '77L', '40D', '86E', '174E', '101D', '162W', '174K', '219W', '224K', '174D', '75I', '75M', '211G', '228Q', '200V', '70N', '39D', '238I', '142T', '215S', '214C', '197L', '54I', '207A', '210W', '138K', '132V', '173I', '214L', '221C', '62V', '68N', '65R', '174R', '188I', '223N', '207K', '210F', '180V', '159K', '141A', '175Y', '219Q', '98G', '191T', '184V', '219D', '139A', '197E', '203K', '172K', '46T', '58N', '194K', '238N', '43N', '60L', '194T', '223Q', '90I', '175H', '32Q', '203G', '101P', '34I', '69A', '39Q', '142V', '83K', '123N', '26M', '16V', '68K', '19S', '187I', '166Q', '70S', '70T', '234I', '101H', '208Y', '200A', '39S', '103R', '21F', '27S', '181C', '75L', '190E', '102Q', '215F', '197P', '67E', '140Q', '174L', '101I', '67Y', '164L', '195K', '103H', '162A', '122I', '73L', '43R', '214I', '102L', '69S', '188H', '210D', '122A', '203R', '208L', '70A', '181G', '203Q', '53D', '211E', '223E', '215E', '31V', '157A', '35I', '43E', '238T', '163I', '221Y', '33G', '75S', '145V', '203A', '215N', '228F', '74I', '70Q', '165I', '60I', '215D', '173S', '199K', '218R', '73M', '44A', '163A', '166T', '48W', '177G', '211A', '32E', '181I', '178M', '165M', '190Q', '190T', '219N', '101E', '100V', '179F', '139M', '194D', '35K', '240S', '135T', '58G', '174H', '173N', '182S', '122D', '36A', '138R', '40F', '11R', '135M', '162T', '39I', '63K', '64R', '182H', '200R', '40K', '240K', '106A', '200M', '219R', '119I', '28K', '68T', '198E', '162F', '107P', '197I', '214M', '197T', '21A', '173D', '13N', '190S', '202T', '207D', '49R', '32R', '224D', '105A', '237E', '135L', '179A', '228R', '211N', '200F', '184L', '90L', '65N', '189I', '73R', '48E', '173R', '142M', '176T', '103G', '173T', '196K', '174G', '238R', '106M', '228H', '31L', '39M', '173L', '88S', '39K', '196E', '177N', '101S', '219E', '70R', '190A', '188C', '46Q', '145E', '201I', '165A', '236L', '170S', '179T', '68R', '67H', '151M', '36D', '208F', '178L', '134G', '228K', '216L', '49Q', '43S', '101R', '28A', '106I', '179L', '123I', '104Q', '89K', '159L', '67N', '162N', '104H', '121C', '163C', '210S', '135R', '127C', '31T', '109V', '64N', '35Q', '69G', '226H', '237N', '39N', '113E', '204A', '43M', '203V', '41L', '200I', '69K', '44D', '23N', '16T', '104R', '205F', '211S', '37L', '204Q', '194R', '144F', '20T', '116Y', '139S', '29K', '170T', '11T', '75A', '66R', '207G', '179E', '123G', '122E', '122R', '43Q', '215C', '195M', '207N', '238Q', '219H', '173E', '139R', '39A', '179I', '230L', '200N', '35M', '87L', '203N', '135K', '145C', '142A', '21I', '102R', '44K', '139K', '62P', '174A', '48V', '196R', '32T', '43A', '207T', '176S', '228S', '70E', '230K', '196Q', '47F', '169H', '69D', '200E', '197K', '174N', '70G', '218E', '194N', '103S', '215L', '223T', '137S', '188L', '69H', '146F', '232H', '49E', '14Q', '169K', '177E', '22R', '64Y', '28G', '178F', '111I', '189L', '163T', '102E', '166I', '167V', '161P', '223R', '45E', '211Q', '35L', '207E', '48T', '132L', '165L', '181S', '36K', '39R', '223A', '200L', '75T', '173A', '28R', '69I', '104N', '102N', '117A', '181F', '115F', '204N', '215Y', '86N', '219G', '14A', '190V', '94L', '197H', '11Q', '135V', '227Y', '40S', '118I', '101N', '122Q', '204K', '50V', '212S', '169D', '207V', '139I', '196D', '27P', '39E', '84V', '217Y', '128A', '79D', '98S', '158S', '203D', '121Y', '138G', '122P', '47V', '96R', '60A', '123S', '215H', '39L', '35T', '101Q', '33V', '171Y')
NRTI_muts<- c('48E', '104N', '37L', '101I', '33G', '101Q', '123I', '227Y', '106A', '21I', '175Y', '173E', '35A', '236R', '167M', '35R', '159K', '39K', '200E', '39N', '139R', '104R', '144F', '157A', '139I', '11T', '60A', '60I', '170T', '45E', '212S', '32T', '14A', '70T', '223Q', '39M', '192N', '170S', '48T', '32E', '190C', '73M', '166I', '67G', '188H', '176S', '158S', '137S', '162F', '216L', '202T', '107A', '13N', '211N', '151M', '39E', '14Q', '215L', '138Q', '44A', '28K', '196R', '67Y', '107S', '28D', '60L', '135K', '102N', '197K', '65R', '219Q', '123E', '219N', '67H', '162H', '101H', '88G', '215V', '16V', '11R', '218E', '111I', '28A', '194D', '179T', '69A', '223N', '102Q', '215F', '200A', '196K', '238R', '162C', '198E', '178F', '69P', '135P', '106I', '105L', '142T', '22N', '174K', '43S', '62V', '90I', '165L', '122P', '36K', '40F', '180V', '211S', '215N', '64R', '223R', '162N', '113E', '11Q', '138G', '77L', '224K', '101D', '132L', '228H', '94L', '39D', '26M', '35K', '98S', '102E', '167V', '103G', '139K', '68T', '205F', '65N', '182S', '103N', '207H', '107P', '169D', '179E', '159V', '221Y', '195K', '46T', '194N', '69I', '135T', '237N', '203N', '69K', '239I', '200M', '207R', '223T', '211I', '101P', '96R', '105A', '135L', '215Y', '237E', '141A', '36D', '181V', '142V', '36A', '179I', '173S', '123Q', '98G', '69N', '173I', '122A', '139S', '210S', '174L', '173D', '211G', '104Q', '219T', '219H', '67E', '70R', '219E', '35T', '204D', '181I', '194T', '39A', '82R', '123S', '136E', '146F', '166Q', '121C', '75L', '173Q', '200V', '163A', '85R', '197I', '207G', '35M', '188I', '122I', '182H', '207E', '68N', '121Y', '70Q', '210F', '104H', '214I', '210M', '62P', '122D', '80I', '194K', '122R', '225H', '101R', '208F', '204N', '228I', '88C', '189L', '211A', '70S', '189I', '181G', '197L', '184I', '195M', '73R', '178L', '230L', '174A', '196E', '227L', '158T', '40S', '101A', '32Q', '40K', '101E', '41I', '163C', '47F', '173T', '106M', '232H', '31L', '100I', '208L', '194R', '188F', '39Q', '86N', '177G', '200F', '162D', '210W', '145C', '207A', '215S', '221C', '240E', '100V', '202V', '195L', '69E', '34I', '165I', '19S', '64H', '211K', '145V', '116Y', '31R', '101N', '214C', '219W', '207D', '174R', '164I', '118C', '108I', '226H', '162W', '166R', '68G', '73L', '64N', '70A', '138K', '211D', '207T', '203D', '203Q', '190E', '54I', '75T', '215H', '103S', '219D', '174G', '16T', '135V', '181F', '39I', '188L', '27P', '29K', '121H', '169K', '174N', '132V', '169H', '176T', '135M', '211E', '75M', '70G', '142A', '214L', '40D', '47V', '43R', '169A', '204A', '27S', '75I', '67N', '145E', '70E', '142K', '74I', '166T', '118I', '204Q', '135R', '89K', '69G', '178M', '190S', '203V', '162T', '238Q', '127C', '43N', '35I', '228R', '69S', '234I', '31V', '123N', '188C', '68K', '13R', '200K', '228N', '175H', '197T', '23N', '177E', '191T', '119I', '224D', '64Y', '212M', '102R', '174E', '43Q', '53D', '88S', '228Q', '142L', '58N', '46Q', '211T', '63K', '228F', '115F', '101S', '203A', '177N', '203G', '103R', '240K', '197E', '49R', '200I', '58G', '238N', '139M', '172K', '20R', '165A', '87L', '215E', '32R', '160L', '174H', '228K', '43E', '44D', '35Q', '238T', '50V', '102L', '48W', '179D', '173A', '48V', '161P', '68R', '163I', '204K', '215I', '128A', '35L', '32N', '28G', '70N', '20T', '179L', '173R', '122Q', '201I', '142M', '215D', '208Y', '140Q', '39L', '230K', '90L', '219R', '122E', '190Q', '211Q', '35E', '67A', '39S', '159L', '203K', '162A', '139A', '179F', '200R', '207K', '197H', '79D', '181C', '109V', '200G', '22R', '66R', '86E', '215C', '102M', '219G', '21A', '28R', '200S', '217Y', '43A', '195T', '223E', '174D', '171Y', '75A', '49E', '223A', '197P', '164L', '33V', '138A', '173N', '184V', '200L', '21F', '162Y', '109I', '74V', '196Q', '178V', '83K', '75S', '196D', '207N', '41L', '69D', '178R', '123G', '44K', '190A', '163T', '69H')
PI_muts<-c('84C', '77G', '54A', '54T', '36A', '70V', '19M', '41E', '11I', '43T', '47V', '98I', '24F', '45N', '82A', '72L', '20R', '72R', '72K', '66V', '63T', '94S', '30N', '69P', '50L', '19V', '83P', '12M', '20L', '77T', '37K', '64V', '38V', '14N', '12R', '66T', '16Q', '34Q', '61N', '63E', '48Q', '99Y', '34S', '41K', '72T', '41P', '89T', '23F', '74P', '63A', '50V', '88G', '36L', '77I', '73S', '95V', '69K', '18I', '22T', '33I', '35A', '90V', '37H', '55T', '67Y', '74S', '73T', '10Y', '71V', '16E', '63D', '20Q', '85M', '18K', '24I', '46V', '34R', '62M', '85V', '21Q', '61R', '61K', '71T', '71F', '64M', '43E', '69R', '51A', '72M', '83S', '12A', '17D', '39Q', '95F', '61G', '84A', '70N', '43Q', '54V', '74K', '79S', '20T', '48S', '12K', '37Y', '14E', '69Q', '82M', '57G', '45V', '71L', '35S', '36T', '22V', '53P', '37T', '48V', '60Q', '37C', '68W', '54M', '53L', '82I', '83H', '67F', '96S', '21D', '72Q', '13V', '37E', '73V', '91N', '19S', '69N', '72E', '33F', '19P', '15R', '59F', '20N', '34A', '35H', '10I', '89V', '70T', '41Q', '64L', '18L', '66L', '60N', '18H', '72N', '47E', '91S', '43I', '74A', '88D', '89F', '65K', '67D', '12I', '19R', '34N', '73A', '41A', '17E', '37V', '41I', '19W', '73I', '83D', '19E', '48I', '74D', '66F', '71M', '67M', '89I', '54S', '11L', '93L', '61D', '53I', '37D', '33M', '88T', '79D', '35D', '70E', '45L', '67S', '72V', '92R', '16A', '39T', '92K', '55R', '77L', '48E', '67W', '10M', '34H', '21K', '10R', '70Q', '79N', '10C', '76V', '36V', '79H', '74Q', '46I', '32V', '47I', '19I', '35N', '34T', '35G', '18R', '18E', '75I', '67G', '15L', '10S', '95L', '65D', '39L', '45Q', '79A', '73D', '76L', '10H', '63S', '67L', '63R', '57K', '84V', '95G', '12D', '15V', '37P', '53Y', '43R', '91A', '10F', '68E', '20M', '14T', '58E', '41N', '20V', '20E', '79Q', '84T', '91K', '36I', '93M', '72F', '88S', '48M', '69Y', '60K', '46L', '89A', '60Y', '32I', '85L', '97I', '80P', '82F', '12Q', '39S', '59H', '12P', '32M', '90M', '48A', '37Q', '16N', '37A', '12S', '19T', '55N', '34K', '54L', '33V', '12V', '71I', '61E', '67E', '15N', '34D', '13L', '77Q', '10V', '45R', '70R', '24M', '14R', '12E', '19Q', '94D', '68D', '62V', '74E', '82T', '92E', '34V', '12N', '72Y', '82L', '37S', '45I', '67H', '89M', '63C', '13A', '43S', '63H', '35Q', '23I', '61H', '82C', '63P', '63V', '38I', '63Q', '82S', '47A', '13M', '60E', '82V', '20I', '73C')

###Non-polymorphic TSMs (derived from Rhee et al. paper)
INIs_tsm <- c('51Y','66A','66I','66K','92Q','92G','95K','97A','118R','121Y','138A','138K','140A','140C','140S','143A','143C','143G','143H','143R','143S','145S','146P','147G','148H','148K','148R','148N','151L','151A','153F','153Y','155H','157Q','163K','232N','263K')
NRTIs_tsm<- c('41L', '43E', '43N', '43Q', '44A', '44D',  '62V', '65R', '67E', '67G', '67N', '69D', '69N', '69S', '70R', '74I', '74V', '75I', '75M', '75T', '77L', '98G', '115F', '116Y', '151M', '184I', '184V', '203K', '208Y', '210W', '215F', '215I', '215V', '215Y', '218E', '219E', '219N', '219Q', '219R', '223Q', '228H', '228R')
NNRTIs_tsm<- c('100I', '101E', '101N', '101P', '103N', '103S', '106A', '106M', '108I', '138Q', '181C', '181I', '181V', '188C', '188H', '188L', '190A', '190E', '190Q', '190S', '221Y', '225H', '227L', '230L', '236L', '238T')
PIs_tsm<- c('10F', '10R', '11I', '20I', '20T', '20V', '23I', '24I', '30N', '32I', '33F', '34Q', '35G', '43T', '46I', '46L', '46V', '47A', '47V',  '48M', '48V', '50L', '50V', '53L', '53Y', '54A', '54L', '54M', '54S', '54T', '54V', '55R', '58E', '66F', '67F', '71I', '73A', '73C', '73S', '73T', '74A', '74P', '74S', '76V', '79A', '82A', '82F', '82S', '82T', '84A', '84C', '84V', '85V', '88D', '88S', '88T', '89V', '90M', '92R', '95F')
#INIs were extracted from the Stanford HIVDB web

###For reproducing the predictions from the report, each dataset should be run the following way:
#INIs
#for(drug in INIs){
#   DRMcv(dataset="INI", drug="RAL", min.muts=2, nfold=5, nrep=10, muts.in = INIs_tsm, add_combinations=FALSE, confusion=FALSE, lars=TRUE)}
#   done
###For getting the predictions with the interactions (LSR-I and LARS-I, it should be run the same with add_combinations=TRUE)

##FOR NRTIs, NNRTIs and PIs, the same applies but for min.muts=10, just change the dataset and drug name