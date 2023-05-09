
## calculate metrics
## library(pROC)
#library(dplyr)

#tmp <- metrics1(label = label,pred_p=pred_p)

#library(purrr)
#
# i <-c(1:100)
#
# a <- list(i)
#
# label1 <- list(label)
# pred_p1 <- list(pred_p)
#
#
# set.seed(123)
# arg1 <- list(label1,pred_p1,c(1:10))
#
# tmp <- arg1 %>% pmap(metrics1)  %>% bind_rows()

#' Title
#'
#' @param label true label
#' @param pred_p predicted prob
#' @param ... XX redundant para

#'
#' @return  mymetrics : dataframe
#' @export
#'
#' @examples


point_metrics <- function(label = NULL,pred_p=NULL,...,i=NULL){


  dat <- data.frame(label = label,
                    pred_p = pred_p)

  ## set positive label; first element is positive,
  # second element is negative
  all_label <- dat %>% group_by(label) %>%
    summarise(pred_mean = mean(pred_p)) %>%
    arrange(-pred_mean) %>%
    pull(label)

  poslabel <- all_label[1]
  neglabel <- all_label[2]

  truelabel <- ifelse(label==poslabel,1,0)


  sample_label <- label
  sample_pred <- pred_p
  pred_label <- ifelse(sample_pred >= 0.5,poslabel,neglabel)

  ## calculate metrics
  ## library(pROC)
  #library(dplyr)
  tmp <- table(sample_label,pred_label ) %>% data.frame()

  auc <- auc(roc(sample_label, sample_pred))[1] %>% round(2)
  tp <- tmp$Freq[tmp[,1]==poslabel & tmp[,2]==poslabel]
  tn <- tmp$Freq[tmp[,1]==neglabel& tmp[,2]==neglabel]
  fp <- tmp$Freq[tmp[,1]==neglabel& tmp[,2]==poslabel]
  fn <- tmp$Freq[tmp[,1]==poslabel& tmp[,2]==neglabel]
  ## sen, recall,
  tpr <- tp/sum(truelabel)
  fnr <- 1- tpr

  tnr <- tn/(nrow(dat)-sum(truelabel))
  fpr <- 1-tnr

  #https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
  prevalence <- (100*mean(truelabel)) %>% round(2)
  acc <- ((tp+tn)/(tp+tn+fp+fn)) %>% round(2)
  sen <- (tp/(tp+fn)) %>% round(2)
  spe <- (tn/(tn+fp)) %>% round(2)
  youden <- round(sen+spe-1,2)
  ppv <- (tp/(tp+fp)) %>% round(2)
  npv <- (tn/(tn+fn)) %>% round(2)
  F1 <- ((2*tp)/(2*tp+fp+fn)) %>% round(2)
  PLR <- round(tpr/fpr,2)# positive likelihood ratio
  NLR <- round(fnr/tnr,2)
  dor <- (PLR/NLR) %>% round(2) #Diagnostic odds ratio
  precision <- (tp/(tp+fp)) %>% round(2)

  for1  <- 1-npv #False omission rate
  fdr <- 1-ppv # False discovery rate
  #MCC <- (sqrt(tpr*tnr*ppv*npv)-sqrt(fnr*fpr*for1*fdr)) %>% round(2)

  mymetrics <- data.frame(
    auc = auc,
    acc = acc,
    sen = sen,
    spe = spe,
    pre = precision,
    youden = youden,
    ppv = ppv,
    npv = npv,
    FOR = for1,
    FDR = fdr,
    F1 = F1,
    PLR = PLR,
    NLR = NLR,
    DOR = dor)

  return(mymetrics)
}




#' Title
#'
#' @param label true label
#' @param pred_p predicted value
#' @param ... redundant
#' @param i iteration
#'
#' @return one bootstrapped point estimates of metrics: dataframe
#' @export
#'
#' @examples


one_metrics <- function(label = NULL,pred_p=NULL,... ,i=NULL){


  dat <- data.frame(label = label,
                    pred_p = pred_p)

  ## set positive label; first element is positive,
  # second element is negative
 all_label <- dat %>% group_by(label) %>%
    summarise(pred_mean = mean(pred_p)) %>%
    arrange(-pred_mean) %>%
    pull(label)

 poslabel <- all_label[1]
 neglabel <- all_label[2]

 truelabel <- ifelse(label==poslabel,1,0)
  ## sample data with replacement
  index <-  sample(length(label),replace = TRUE)
  sample_label <- label[index]
  sample_pred <- pred_p[index]
  pred_label <- ifelse(sample_pred >= 0.5,poslabel,neglabel)


  tmp <- table(sample_label,pred_label ) %>% data.frame()

  auc <- auc(roc(sample_label, sample_pred))[1] %>% round(2)
  tp <- tmp$Freq[tmp[,1]==poslabel & tmp[,2]==poslabel]
  tn <- tmp$Freq[tmp[,1]==neglabel& tmp[,2]==neglabel]
  fp <- tmp$Freq[tmp[,1]==neglabel& tmp[,2]==poslabel]
  fn <- tmp$Freq[tmp[,1]==poslabel& tmp[,2]==neglabel]
  ## sen, recall,
  tpr <- tp/sum(truelabel)
  fnr <- 1- tpr

  tnr <- tn/(nrow(dat)-sum(truelabel))
  fpr <- 1-tnr

  #https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
  prevalence <- (100*mean(truelabel)) %>% round(2)
  acc <- ((tp+tn)/(tp+tn+fp+fn)) %>% round(2)
  sen <- (tp/(tp+fn)) %>% round(2)
  spe <- (tn/(tn+fp)) %>% round(2)
  youden <- round(sen+spe-1,2)
  ppv <- (tp/(tp+fp)) %>% round(2)
  npv <- (tn/(tn+fn)) %>% round(2)
  F1 <- ((2*tp)/(2*tp+fp+fn)) %>% round(2)
  PLR <- round(tpr/fpr,2)# positive likelihood ratio
  NLR <- round(fnr/tnr,2)
  dor <- (PLR/NLR) %>% round(2) #Diagnostic odds ratio
  precision <- (tp/(tp+fp)) %>% round(2)

  for1  <- 1-npv #False omission rate
  fdr <- 1-ppv # False discovery rate
  #MCC <- (sqrt(tpr*tnr*ppv*npv)-sqrt(fnr*fpr*for1*fdr)) %>% round(2)

  mymetrics <- data.frame(
                          auc = auc,
                          acc = acc,
                          sen = sen,
                          spe = spe,
                          pre = precision,
                          youden = youden,
                          ppv = ppv,
                          npv = npv,
                          FOR = for1,
                          FDR = fdr,
                          F1 = F1,
                          PLR = PLR,
                          NLR = NLR,
                          DOR = dor)

  return(mymetrics)
}





#' Title
#'
#' @param label true label
#' @param pred_p predicted prob
#' @param iteration # of bootstrap iteration, default=1000
#' @param ... redundant para if any
#'
#' @return point and 95% CI dataframe
#' @export
#'
#' @examples
metrics1 <- function(label = NULL,pred_p = NULL,iteration=1000,...){
  set.seed(1017)
  label1 <- list(label)
  pred_p1 <- list(pred_p)
  arg_list <- list(label1,pred_p1,seq_len(iteration))
  tmp_res <- arg_list %>% pmap(one_metrics)%>% bind_rows()
  quantile025 <- tmp_res[,c(1:13)] %>%
    map(~quantile(x = .,probs = 0.025)) %>% unlist()
  quantile975 <-  tmp_res[,c(1:13)] %>%
    map(~quantile(x = .,probs = 0.975)) %>% unlist()

  point_est <- point_metrics(label = label,pred_p = pred_p)

  res <- paste0(point_est,"(",paste(quantile025 %>% round(2),quantile975 %>% round(2),sep = ','),")") %>% data.frame() %>%
    t() %>% data.frame()
  colnames(res) <- colnames(point_est)

  return(res)




}

# #tmp <- metrics1(label = label,pred_p= pred_p,iteration=100)
#
# et_test_proba <- read_excel("E:/2022/AECOPD/from wang/modelling/output/et_test_proba.xlsx")
#
# label <- et_test_proba %>% pull(5)
#
# pred_p <-  et_test_proba %>% pull(3)
#
# roc(label,pred_p)
#
#
# tmp <- metrics1(label,pred_p,iteration = 100)



