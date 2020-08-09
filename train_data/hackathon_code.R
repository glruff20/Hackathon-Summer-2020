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
plot(fit_cd4)

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

#when put this data in, replaced my_predictions2.csv

#issue - still havent included nasal_microbio data into regression model
#another issue - probably want to do a regression with all variables together, that way more important markers will be represented that way in regression?
#could always just add more features to this method too, see if it works any better

save(prediction, cd4_predict, cd8_predict, cd19_predict, nasal_gene_predict, nasal_predict_median, cd4_features, cd8_features, cd19_features, nasal_gene_features, fit_cd19, fit_cd4, fit_cd8, fit_nasal_gene, file="hackathon.rda")

#combine all data into a table, colnames prefix is cell type
head(severe_nasal_micro[,1:5])
colnames(severe_cd19)=c("severity_score", paste("cd19", colnames(severe_cd19[,-c(1)]), sep="_"))
colnames(severe_cd4)=c("severity_score", paste("cd4", colnames(severe_cd4[,-c(1)]), sep="_"))
colnames(severe_cd8)=c("severity_score", paste("cd8", colnames(severe_cd8[,-c(1)]), sep="_"))
colnames(severe_nasal_gene)=c("severity_score", paste("nasal", colnames(severe_nasal_gene[,-c(1)]), sep="_"))
colnames(severe_nasal_micro)=c("severity_score", paste("nasal", colnames(severe_nasal_micro[,-c(1)]), sep="_"))

test=merge(severe_cd19, severe_cd4[,-c(1)], by=0, all.x=T, all.y=T)
head(test2[,1:5])
rownames(test)=test$Row.names
test2=merge(test[,-c(1)], severe_cd8[,-c(1)], by=0, all.x=T, all.y=T)
rownames(test2)=test2$Row.names
test3=merge(test2[,-c(1)], severe_nasal_gene[,-c(1)], by=0, all.x=T, all.y=T)
rownames(test3)=test3$Row.names
test4=merge(test3[,-c(1)], severe_nasal_micro[,-c(1)], by=0, all.x=T, all.y=T)
rownames(test4)=test4$Row.names
test4=test4[,-c(1)]
head(test4[,1:5])
test4[1,]

merged_data=test4 #severity is first column, rest of columns is rest of data

#issue - only have severity scores for patients w/ cd19 data
severe_cd19[,1]
test=merge(severe_cd19[,1:2], severe_cd4[,1:2], by=0, all.x=T, all.y=T)
rownames(test)=test$Row.names
test=test[,-c(1,3,5)]
test
test2=merge(test, severe_cd8[,1:2], by=0, all.x=T, all.y=T)
rownames(test2)=test2$Row.names
test2=test2[,-c(1,5)]
test2
test3=merge(test2, severe_nasal_gene[,1:2], by=0, all.x=T, all.y=T)
rownames(test3)=test3$Row.names
test3
test3=test3[,-c(1,6)]
test4=merge(test3, severe_nasal_micro[,1:2], by=0, all.x=T, all.y=T)
rownames(test4)=test4$Row.names
test4=test4[,-c(1,7)]
test4

merged_severity_score=test4
test4=merge(merged_severity_score, rowMeans(merged_severity_score, na.rm=T), by=0, all.x=T, all.y=T)
rownames(test4)=test4$Row.names
test4=test4[,-c(1)]
test4
merged_severity_score=test4

head(merged_data[,1:5])
test=cbind.data.frame(merged_severity_score[,6], merged_data)
head(test[,1:5])
test=test[,-c(2)]
colnames(test)=c("severity_score", colnames(test[,-c(1)]))

merged_data=test

#calculate correlations for each gene in gene expr tables, subset out top genes to use in regression model
cor_merged=cor(merged_data, use="pairwise.complete.obs")
summary(abs(cor_merged[,1]))
quantile(abs(cor_merged[,1]), 0.999, na.rm=T)

sig_cor_merged=subset(cor_merged, abs(cor_merged[1,]) > 0.55)
dim(sig_cor_merged)
rownames(sig_cor_merged)
sig_cor_merged[,1]


merged_features=subset(merged_data, select=rownames(sig_cor_merged))
head(merged_features)

fit_merged=lm(merged_features$severity_score ~., data=merged_features)
summary(fit_merged)
plot(fit_cd4)

#problem - lots of NA's, would only be able to predict severity for patients w/ data for cd4 or cd8
#select 10 features from each cell type to use, combine those into table

rownames(sig_cor_severe_cd19)=c("severity_score", paste("cd19", rownames(sig_cor_severe_cd19[-c(1),]), sep="_"))
rownames(sig_cor_severe_cd4)=c("severity_score", paste("cd4", rownames(sig_cor_severe_cd4[-c(1),]), sep="_"))
rownames(sig_cor_severe_cd8)=c("severity_score", paste("cd8", rownames(sig_cor_severe_cd8[-c(1),]), sep="_"))
rownames(sig_cor_severe_nasal_gene)=c("severity_score", paste("nasal", rownames(sig_cor_severe_nasal_gene[-c(1),]), sep="_"))


rownames(sig_cor_severe_cd8)
test_features=subset(merged_data, select=c(rownames(sig_cor_severe_cd19[1:28,]), rownames(sig_cor_severe_cd4[-c(1),]), rownames(sig_cor_severe_cd8[-c(1),]), rownames(sig_cor_severe_nasal_gene[-c(1),])))
head(test_features)

fit_merged=lm(test_features$severity_score ~., data=test_features)
summary(fit_merged)


#replace blank values with median score that gene
library(tidyr)
col_means <- lapply(test_features, median, na.rm = TRUE)
test1 <- replace_na(test_features, col_means)
head(test1)
test_features_nona=test1

fit_merged=lm(test_features_nona$severity_score ~., data=test_features_nona)
summary(fit_merged)

merged_cor=cor(test_features_nona, use="pairwise.complete.obs")
summary(abs(merged_cor[,1]))
quantile(abs(merged_cor[,1]), 0.85)


sig_cor_merged=subset(merged_cor, abs(merged_cor[1,]) > 0.405)  #has the highest adjusted r-squared value
dim(sig_cor_merged)
rownames(sig_cor_merged)
sig_cor_merged[,1]


test=subset(test_features_nona, select=rownames(sig_cor_merged))
head(test)

fit_merged=lm(severity_score ~., data=test)
summary(fit_merged)

fit_test=lm(severity_score ~., data=sig_features)
summary(fit_merged)
plot(fit_merged)
plot(fit_merged$residuals)

sig_features=test

library(DAAG)

cross_v=cv.lm(test, fit_merged, m=3)
summary(cross_v)
head(cross_v)
cross_v[,c(1,62:63)]

pr=residuals(fit_merged)/(1-lm.influence(fit_merged)$hat)
press=sum(pr^2)

my.anova=anova(fit_merged)
tss=sum(my.anova$`Sum Sq`)
pred.r.square=1-(press/(tss))
pred.r.square
#cross validation looks pretty good, need to add in nasal_microbio data
#for microbio data, just do 1's and 0's at first

#turn nasal micro data set into binary, 1 or 0
head(severe_nasal_micro[,1:5])
test=severe_nasal_micro
test[test>0.0000001] = 1
head(test[,1:5])
nasal_micro_bi=cbind.data.frame(severe_nasal_micro[,1], test[,-c(1)])
colnames(nasal_micro_bi)=c("severity_score", colnames(nasal_micro_bi[,-c(1)]))
head(nasal_micro_bi[,1:10])

#merge nasal micro into big dataset from before

head(merged_data[,1:5])
test4=merge(test3[,-c(1)], nasal_micro_bi[,-c(1)], by=0, all.x=T, all.y=T)
rownames(test4)=test4$Row.names
test4=test4[,-c(1)]

test=cbind.data.frame(merged_severity_score[,6], test4)
head(test[,1:5])
test=test[,-c(2)]
colnames(test)=c("severity_score", colnames(test[,-c(1)]))

merged_bi_nasal=test
tail(merged_bi_nasal[,24650:24658])

#calculate correlations for nasal_micro dataset
#take out top 20 or so genes

cor_nasal_micro=cor(nasal_micro_bi, use="pairwise.complete.obs")
summary(abs(cor_nasal_micro[,1]))
quantile(abs(cor_nasal_micro[,1]), 0.87, na.rm=T)
sig_cor_nasal_micro=subset(cor_nasal_micro, abs(cor_nasal_micro[1,]) > 0.19)
dim(sig_cor_nasal_micro)
rownames(sig_cor_nasal_micro)
sig_cor_severe_nasal_gene[,1]

test_features1=subset(merged_bi_nasal, select=c(rownames(sig_cor_severe_cd19[1:28,]), rownames(sig_cor_severe_cd4[-c(1),]), rownames(sig_cor_severe_cd8[-c(1),]), rownames(sig_cor_severe_nasal_gene[-c(1),]), rownames(sig_cor_nasal_micro[-c(1),])))
head(test_features1)

head(test_features1[,1:135])
col_median=lapply(test_features1[,1:135], median, na.rm=T)
test6=replace_na(test_features1[,1:135], col_median)
###may need to add replacements to NA's in the microbio columns
head(test6)
test7=cbind.data.frame(test6, test_features1[,136:152])
head(test7)
test_features_nona1=test7


merged_cor1=cor(test_features_nona1, use="pairwise.complete.obs")
summary(abs(merged_cor1[,1]))
quantile(abs(merged_cor1[,1]), 0.85)

sig_cor_merged1=subset(merged_cor1, abs(merged_cor1[1,]) > 0.355)  #has the highest adjusted r-squared value
dim(sig_cor_merged1)
rownames(sig_cor_merged1)
sig_cor_merged[,1]


test=subset(test_features_nona1, select=rownames(sig_cor_merged1))
head(test)

fit_merged1=lm(severity_score ~., data=test)
summary(fit_merged1)

###the top adjusted p value had no features for nasal microbiome data
###there is one patient in the test data taht only has nasal microbio data
### for them maybe make a separate regression model for only microbio data

##try using stepwise regression to pick best regression model from dataframes

test_features_nona1 #150 features correlated with severity score
library(MASS)
full_model=lm(severity_score ~.,data=test)
step.model=stepAIC(full_model, direction="both", trace=F)
summary(step.model)

pred.r.squared <- pred_r_squared(step.model)
pred.r.squared

library(caret)
train.control=trainControl(method="cv", number=10)
step.model=train(severity_score~., data=test_features_nona, method="leapSeq", trControl=train.control)
step.model$results
step.model$bestTune
summary(step.model$finalModel)
coef(step.model$finalModel, 4)
test_model=lm(severity_score~cd4_ENSG00000130638.15 + cd8_ENSG00000115053.15 + nasal_ENSG00000198805.11 + cd8_ENSG00000075415.12, data=test_features_nona)
summary(test_model)

pred.r.squared <- pred_r_squared(test_model)
pred.r.squared

summary(fit_merged)
