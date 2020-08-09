#running the model on test data

#download test data tables
test_cd19=read.table("cd19_gene_expr_test.txt", check.names=F)
test_cd4=read.table("cd4_gene_expr_test.txt", check.names=F)
test_cd8=read.table("cd8_gene_expr_test.txt", check.names=F)
test_nasal_gene=read.table("nasal_gene_expr_test.txt", check.names=F)

#subset out useful columns for each data file
test_cd19=t(test_cd19)
test_cd4=t(test_cd4)
test_cd8=t(test_cd8)
test_nasal_gene=t(test_nasal_gene)

head(colnames(test_cd19))
rownames(sig_cor_severe_cd19)
test_cd19_features=subset(test_cd19, select=rownames(sig_cor_severe_cd19))  #subscript out of bounds, delete column for severity score now
test_cd19_features=subset(test_cd19, select=rownames(sig_cor_severe_cd19[2:28,]))
head(test_cd19_features)

rownames(sig_cor_severe_cd8)
test_cd8_features=subset(test_cd8, select=rownames(sig_cor_severe_cd8[2:15,]))
head(test_cd8_features)

rownames(sig_cor_severe_cd4)
test_cd4_features=subset(test_cd4, select=rownames(sig_cor_severe_cd4[2:48,]))
head(test_cd4_features)

rownames(sig_cor_severe_nasal_gene)
test_nasal_gene_features=subset(test_nasal_gene, select=rownames(sig_cor_severe_nasal_gene[2:47,]))
head(test_nasal_gene_features)

#use linear regression models to make predictions on test data
cd19_predict=predict(fit_cd19, data.frame(test_cd19_features)) #need to coerce new data table into dataframe - checked to make sure no info lost
cd4_predict=predict(fit_cd4, data.frame(test_cd4_features))
cd8_predict=predict(fit_cd8, data.frame(test_cd8_features))
nasal_gene_predict=predict(fit_nasal_gene, data.frame(test_nasal_gene_features))

#combine predictions into table
write.table(cd4_predict, file="cd4_predictions.txt")
write.table(cd8_predict, file="cd8_predictions.txt")
write.table(cd19_predict, file="cd19_predictions.txt")
write.table(nasal_gene_predict, file="nasal_gene_predictions.txt")

cd4_predict=read.table("cd4_predictions.txt")
cd8_predict=read.table("cd8_predictions.txt")
cd19_predict=read.table("cd19_predictions.txt")
nasal_gene_predict=read.table("nasal_gene_predictions.txt")


cd4_cd8=merge(cd4_predict, cd8_predict, by=0, all.x=T, all.y=T)
rownames(cd4_cd8)=cd4_cd8$Row.names
cd4_cd8=cd4_cd8[,-c(1)]
cd4_8_19=merge(cd4_cd8, cd19_predict, by=0, all.x=T, all.y=T)
head(cd4_8_19)
rownames(cd4_8_19)=cd4_8_19$Row.names
cd4_8_19=cd4_8_19[,-c(1)]
merged_predict=merge(cd4_8_19, nasal_gene_predict, by=0, all.x=T, all.y=T)
rownames(merged_predict)=merged_predict$Row.names
merged_predict=merged_predict[,-c(1)]
colnames(merged_predict)=c("cd4", "cd8", "cd19", "nasal_gene")
predict_average=rowMeans(merged_predict, na.rm=T)
test=cbind.data.frame(merged_predict, predict_average)
merged_predict=test
head(merged_predict)

#match numbers up with their order for prediction.csv file

prediction=read.table("prediction.csv", sep=",")
head(prediction)
prediction=prediction[-c(1),]
rownames(prediction)=prediction$V2
test=merge(prediction, merged_predict, by=0, all.x=T, all.y=T, sort=F)
test

test2=cbind.data.frame(test$predict_average, test$Row.names)
head(test2)

write.table(test2, file="my_predictions.csv", sep=",") #first draft, problems below


#issues with current predictions - nasal gene large variance
#3 patients dont have any predictions (no data) - at bottom of file right now
#figure out how to make nasal gene predictions better - just use median?
#how to use nasal microbio dataset??

#first attempt to fix - just use median for nasal gene predictions (only for patients that have nasal gene data)
nasal_predict_median=cbind.data.frame(nasal_gene_predict, rep(median(test$nasal_gene, na.rm=T), 42))

merged_predict2=merge(cd4_8_19, nasal_predict_median, by=0, all.x=T, all.y=T)
rownames(merged_predict2)=merged_predict2$Row.names
head(merged_predict2)
merged_predict2=merged_predict2[,-c(1,5)]
colnames(merged_predict2)=c("cd4", "cd8", "cd19", "nasal_gene_median")
predict_average2=rowMeans(merged_predict2, na.rm=T)
merged_predict2=cbind.data.frame(merged_predict2, predict_average2)

test=merge(prediction, merged_predict2, by=0, all.x=T, all.y=T, sort=F)

test2=cbind.data.frame(test$predict_average, test$Row.names)
head(test2)

write.table(test2, file="my_predictions2.csv", sep=",") #use median of nasal_gene table, still missing 3 patients severity scores

#fix missing patients, determine what patients in training data dont have data, calculate average of those, plug in that average for patients in test data

not_nasal_gene=subset(severity_score, !(rownames(severity_score) %in% rownames(nasal_gene)))
not_nasal_cd4=subset(not_nasal_gene, !(rownames(not_nasal_gene) %in% rownames(cd4_gene)))
severity_no_data=not_nasal_cd4 #2 very low severity scores and 1 high severity score - doesnt seem very predictive...
mean(as.numeric(paste(severity_no_data$V1))) # mean is 3.1768 for patients with no data
#put this data in by hand, since still need to add nasal_microbio data for one of the test patients

#when put this data in, replaced my_predictions2.csv


#running new model - fit_merged
test_cd19=read.table("cd19_gene_expr_test.txt", check.names=F)
test_cd4=read.table("cd4_gene_expr_test.txt", check.names=F)
test_cd8=read.table("cd8_gene_expr_test.txt", check.names=F)
test_nasal_gene=read.table("nasal_gene_expr_test.txt", check.names=F)
test_nasal_micro=read.table("nasal_microbiome_test.txt", check.names=F)

test_cd19=t(test_cd19)
test_cd4=t(test_cd4)
test_cd8=t(test_cd8)
test_nasal_gene=t(test_nasal_gene)
test_nasal_micro=t(test_nasal_micro)

colnames(test_cd19)=paste("cd19", colnames(test_cd19), sep="_")
colnames(test_cd4)=paste("cd4", colnames(test_cd4), sep="_")
colnames(test_cd8)=paste("cd8", colnames(test_cd8), sep="_")
colnames(test_nasal_gene)=paste("nasal", colnames(test_nasal_gene), sep="_")
colnames(test_nasal_micro)=paste("nasal", colnames(test_nasal_micro), sep="_")

#merge cell type datas together into big table

test=merge(test_cd19, test_cd4, by=0, all.x=T, all.y=T)
head(test[,1:5])
rownames(test)=test$Row.names
test2=merge(test[,-c(1)], test_cd8, by=0, all.x=T, all.y=T)
rownames(test2)=test2$Row.names
test3=merge(test2[,-c(1)], test_nasal_gene, by=0, all.x=T, all.y=T)
rownames(test3)=test3$Row.names
test4=merge(test3[,-c(1)], test_nasal_micro, by=0, all.x=T, all.y=T)
rownames(test4)=test4$Row.names
test4=test4[,-c(1)]
head(test4[,1:5])

merged_test_data=test4 #severity is first column, rest of columns is rest of data

#subset out only features used in regression model

merged_test_features=subset(merged_test_data, select=rownames(sig_cor_merged[-c(1),]))
head(merged_test_features)

#get rid of NA values in this dataset, replace with median cell values
library(tidyr)
test_col_medians <- lapply(merged_test_features, median, na.rm = TRUE)
test1 <- replace_na(merged_test_features, test_col_medians)
head(test1)
merged_test_features_nona=test1

#predict severity score with regression model
#fit_merged has predicted r-squared value of 0.44
summary(fit_merged)
test_predictions=predict(fit_merged, data.frame(merged_test_features_nona))
test_predictions
test5=cbind.data.frame(test_predictions, test_predictions)
test5
test_predictions=test5


#create functions to calculate predicted r-squared value for linear regressions
pred_r_squared <- function(linear.model) {
  lm.anova <- anova(linear.model)
  tss <- sum(lm.anova$"Sum Sq")
  # predictive R^2
  pred.r.squared <- 1 - PRESS(linear.model)/(tss)
  return(pred.r.squared)
}

PRESS <- function(linear.model) {
  pr <- residuals(linear.model)/(1 - lm.influence(linear.model)$hat)
  PRESS <- sum(pr^2)
  return(PRESS)
}

pred.r.squared <- pred_r_squared(fit_merged)
pred.r.squared

###match predictions up to prediction file
prediction=read.table("prediction.csv", sep=",")
head(prediction)
prediction=prediction[-c(1),]
rownames(prediction)=prediction$V2
test=merge(prediction, test_predictions, by=0, all.x=T, all.y=T, sort=F)
test

test2=cbind.data.frame(test$test_predictions, test$Row.names)
head(test2)

write.table(test2, file="my_predictions3.csv", sep=",") #first draft, problems below


#with step.model - from stepAIC - predicted r squared good, predictions awful...?

predict(step.model, data.frame(test5))

step.model$model
test5=subset(merged_test_data, select=colnames(test[,-c(1)]))
head(test5)

#get rid of NA values in this dataset, replace with median cell values
library(tidyr)
test_col_medians <- lapply(test5, median, na.rm = TRUE)
test5 <- replace_na(test5, test_col_medians)
head(test5)
merged_test_features_nona=test1

        