#hackathon-urmc-2020

#download datasets
severity_score=read.table("severity_score_train.txt")
nasal_micro=read.table("nasal_microbiome_train.txt", check.names=F)
nasal_gene=read.table("nasal_gene_expr_train.txt", check.names=F)
cd8_gene=read.table("cd8_gene_expr_train.txt", check.names=F)
cd4_gene=read.table("cd4_gene_expr_train.txt", check.names=F)
cd19_gene=read.table("cd19_gene_expr_train.txt", check.names=F)

#create tables of severity score as column 1, genes as rest of columns
#rows are all the patient ids

severity_score=severity_score[-c(1),]
rownames(severity_score)=severity_score[,2]

t_severe=t(severity_score)
nasal_gene=t(nasal_gene)
test=merge(severity_score, nasal_gene, by=0)
head(test[,1:5])
rownames(test)=test[,1]
test=test[,-c(1)]
test=test[,-c(2)]
colnames(test)=c("severity_score", colnames(test[,2:6845]))
severe_nasal_gene=test

nasal_micro=t(nasal_micro)
severe_nasal_micro=merge(severity_score, nasal_micro, by=0)
head(severe_nasal_micro[,1:5])
rownames(severe_nasal_micro)=severe_nasal_micro[,1]
severe_nasal_micro=severe_nasal_micro[,-c(1,3)]
colnames(severe_nasal_micro)=c("severity_score", colnames(severe_nasal_micro[,2:149]))

cd8_gene=t(cd8_gene)
severe_cd8=merge(severity_score, cd8_gene, by=0)
head(severe_cd8[,1:5])
rownames(severe_cd8)=severe_cd8[,1]
severe_cd8=severe_cd8[,-c(1,3)]
colnames(severe_cd8)=c("severity_score", colnames(severe_cd8[,2:5828]))

cd4_gene=t(cd4_gene)
severe_cd4=merge(severity_score, cd4_gene, by=0)
head(severe_cd4[,1:5])
rownames(severe_cd4)=severe_cd4[,1]
severe_cd4=severe_cd4[,-c(1,3)]
colnames(severe_cd4)=c("severity_score", colnames(severe_cd4[,2:5813]))

cd19_gene=t(cd19_gene)
severe_cd19=merge(severity_score, cd19_gene, by=0)
head(severe_cd19[,1:5])
rownames(severe_cd19)=severe_cd19[,1]
severe_cd19=severe_cd19[,-c(1,3)]
colnames(severe_cd19)=c("severity_score", colnames(severe_cd19[,2:5828]))

#calculate correlation scores for each gene vs severity score; take top portion

#make severity score column numeric for all data tables
severe_cd19[,1]=as.numeric(paste(severe_cd19[,1]))
severe_cd4[,1]=as.numeric(paste(severe_cd4[,1]))
severe_cd8[,1]=as.numeric(paste(severe_cd8[,1]))
severe_nasal_gene[,1]=as.numeric(paste(severe_nasal_gene[,1]))
severe_nasal_micro[,1]=as.numeric(paste(severe_nasal_micro[,1]))

#calculate correlations for each gene in gene expr tables, subset out top genes to use in regression model
cor_severe_cd19=cor(severe_cd19, use="pairwise.complete.obs")
summary(cor_severe_cd19[,1])
sig_cor_severe_cd19=subset(cor_severe_cd19, cor_severe_cd19[1,] < -0.43 | cor_severe_cd19[1,] > 0.43)
dim(sig_cor_severe_cd19)
rownames(sig_cor_severe_cd19)
sig_cor_severe_cd19[,1]

cor_severe_cd8=cor(severe_cd8, use="pairwise.complete.obs")
summary(cor_severe_cd8[,1])
sig_cor_severe_cd8=subset(cor_severe_cd8, cor_severe_cd8[1,] < -0.55 | cor_severe_cd8[1,] > 0.6)
dim(sig_cor_severe_cd8)
rownames(sig_cor_severe_cd8)
sig_cor_severe_cd8[,1]

cor_severe_cd4=cor(severe_cd4, use="pairwise.complete.obs")
summary(cor_severe_cd4[,1])
sig_cor_severe_cd4=subset(cor_severe_cd4, cor_severe_cd4[1,] < -0.35 | cor_severe_cd4[1,] > 0.4)
dim(sig_cor_severe_cd4)
rownames(sig_cor_severe_cd4)
sig_cor_severe_cd4[,1]

cor_severe_nasal_gene=cor(severe_nasal_gene, use="pairwise.complete.obs")
summary(cor_severe_nasal_gene[,1])
sig_cor_severe_nasal_gene=subset(cor_severe_nasal_gene, cor_severe_nasal_gene[1,] < -0.40 | cor_severe_nasal_gene[1,] > 0.40)
dim(sig_cor_severe_nasal_gene)
rownames(sig_cor_severe_nasal_gene)
sig_cor_severe_nasal_gene[,1]



#figure out how to correlate nasal micro data??
cor_severe_nasal_micro=cor(severe_nasal_micro, use="pairwise.complete.obs")
summary(cor_severe_nasal_gene[,1])
sig_cor_severe_nasal_gene=subset(cor_severe_nasal_gene, cor_severe_nasal_gene[1,] < -0.43 | cor_severe_nasal_gene[1,] > 0.45)
dim(sig_cor_severe_nasal_gene)
rownames(sig_cor_severe_nasal_gene)
sig_cor_severe_nasal_gene[,1]


#now that have features chosen, subset out tables to have severity as col1
#all significant features as the rest of the columns

cd4_features=subset(severe_cd4, select=rownames(sig_cor_severe_cd4))
cd8_features=subset(severe_cd8, select=rownames(sig_cor_severe_cd8))
sig_cor_severe_cd19=sig_cor_severe_cd19[1:28,] #29th row was NA, not sure why, deleted it
cd19_features=subset(severe_cd19, select=rownames(sig_cor_severe_cd19))
nasal_gene_features=subset(severe_nasal_gene, select=rownames(sig_cor_severe_nasal_gene))


#create a linear regression model from those data tables for each dataset
library(DAAG) #use for cross validation testing

fit_cd4=lm(cd4_features$severity_score ~., data=cd4_features)
summary(fit_cd4)

fit_cd8=lm(cd8_features$severity_score ~., data=cd8_features)
summary(fit_cd8)

fit_cd19=lm(cd19_features$severity_score ~., data=cd19_features)
summary(fit_cd19)

fit_nasal_gene=lm(nasal_gene_features$severity_score ~., data=nasal_gene_features)
summary(fit_nasal_gene)

#now have regression models for each data table; use these to predict severity scores for each patient in test data

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