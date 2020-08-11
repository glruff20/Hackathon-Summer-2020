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

#figure out how to correlate nasal micro data?? - did later in code
#cor_severe_nasal_micro=cor(severe_nasal_micro, use="pairwise.complete.obs")
#summary(cor_severe_nasal_gene[,1])
#sig_cor_severe_nasal_gene=subset(cor_severe_nasal_gene, cor_severe_nasal_gene[1,] < -0.43 | cor_severe_nasal_gene[1,] > 0.45)
#dim(sig_cor_severe_nasal_gene)
#rownames(sig_cor_severe_nasal_gene)
#sig_cor_severe_nasal_gene[,1]


#now that have features chosen, subset out tables to have severity as col1
#all significant features as the rest of the columns

cd4_features=subset(severe_cd4, select=rownames(sig_cor_severe_cd4))
cd8_features=subset(severe_cd8, select=rownames(sig_cor_severe_cd8))
sig_cor_severe_cd19=sig_cor_severe_cd19[1:28,] #29th row was NA, not sure why, deleted it
cd19_features=subset(severe_cd19, select=rownames(sig_cor_severe_cd19))
nasal_gene_features=subset(severe_nasal_gene, select=rownames(sig_cor_severe_nasal_gene))


#create a linear regression model from those data tables for each dataset
fit_cd4=lm(cd4_features$severity_score ~., data=cd4_features)
summary(fit_cd4)
plot(fit_cd4)


colnames(cd8_features)=c("severity_score", paste("cd8", colnames(cd8_features[,-c(1)]), sep="_"))
fit_cd8=lm(cd8_features$severity_score ~., data=cd8_features)
summary(fit_cd8)

fit_cd19=lm(cd19_features$severity_score ~., data=cd19_features)
summary(fit_cd19)

fit_nasal_gene=lm(nasal_gene_features$severity_score ~., data=nasal_gene_features)
summary(fit_nasal_gene)

#now have regression models for each data table; use these to predict severity scores for each patient in test data

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

#see if patients show differences btwn each other if missing data
no_cd8=subset(merged_data, is.na(merged_data$cd8_ENSG00000000419.12))
histogram(no_cd8[,1]) #normal median
no_cd4=subset(merged_data, is.na(merged_data$cd4_ENSG00000000419.12))
histogram(no_cd4[,1]) #normal median
no_cd19=subset(merged_data, is.na(merged_data$cd19_ENSG00000000419.12))
summary(no_cd19[,1]) #normal median
no_nasal_gene=subset(merged_data, is.na(merged_data$nasal_ENSG00000000003.14))
histogram(no_nasal_gene[,1], breaks=8) #median about 6
no_nasal_micro=subset(merged_data, is.na(merged_data$`nasal_g__Actinomyces,s__`))
histogram(no_nasal_micro[,1], breaks=8) #median about 6

no_nasal=subset(merged_data, is.na(merged_data$nasal_ENSG00000000003.14) & is.na(merged_data$`nasal_g__Actinomyces,s__`))
histogram(no_nasal[,1], breaks=7) #median about 6

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
plot(severity_score~., data=test)

fit_merged=lm(severity_score ~., data=test)
summary(fit_merged)

fit_merged=lm(severity_score ~., data=sig_features)
summary(fit_merged)
plot(fit_merged)
plot(fit_merged$residuals)

pred.r.squared=pred_r_squared(fit_merged)
pred.r.squared


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

#other things to try: see if regression models from only one cell type are better than combined?

save(merged_data, merged_bi_nasal, merged_test_data, file="merged_data_tables.rda")


#separate patients based on severity
summary(merged_data$severity_score)
low_data=subset(merged_data, merged_data$severity_score<2.5)
med_data=subset(merged_data, merged_data$severity_score>2.5 & merged_data$severity_score<6.5)
high_data=subset(merged_data, merged_data$severity_score>6.5)
head(low_data[,1:5])
averages_data=colMeans(merged_data, na.rm=T)
head(averages_data)

logfc_low=log2(colMeans(low_data, na.rm=T)/averages_data)
head(logfc_low)
summary(logfc_low)
quantile(abs(logfc_low), 0.995, na.rm=T)
change_low=subset(logfc_low, abs(logfc_low)>1.98 & logfc_low !=-Inf)
change_low

logfc_med=log2(colMeans(med_data, na.rm=T)/averages_data)
summary(logfc_med)

logfc_high=log2(colMeans(high_data, na.rm=T)/averages_data)
summary(logfc_high)
quantile(abs(logfc_high), 0.995, na.rm=T)
change_high=subset(logfc_high, abs(logfc_high)>2.12)
head(change_high)
change_high=subset(change_high, change_high !=-Inf)
change_high

change_features=subset(merged_data, select=c("severity_score", names(change_high), names(change_low)))
head(change_features[,1:5])

col_median=lapply(change_features, median, na.rm=T)
test6=replace_na(change_features, col_median)

change_cor=cor(test6)
head(change_cor[,1:5])
summary(change_cor[,1])
sig_change_cor=subset(change_cor, abs(change_cor[,1]) > 0.28)
sig_change_cor=subset(sig_change_cor, rownames(sig_change_cor)!="severity_score.1")
sig_change_features=subset(test6, select=rownames(sig_change_cor))
change_model=lm(severity_score~., data=sig_change_features)
summary(change_model)

head(sig_change_features)

col_means <- lapply(merged_data, median, na.rm = TRUE)
test1 <- replace_na(merged_data, col_means)
head(test1)
test_features_nona=test1

test2=subset(test1, select=c(rownames(sig_cor_merged), "cd19_ENSG00000252614.1", "nasal_ENSG00000163958.13"))
head(test2)
plot(severity_score~., data=test2)

fit_test=lm(severity_score ~., data=test2)
summary(fit_test)

pred.r.squared <- pred_r_squared(fit_test)
pred.r.squared

##adding values from greatest change btwn high and low severity didn't increase predicted r-square

save(fit_merged, file="fit_merged.rda")


#normalize merged data by z-score; make new regression model from that

head(merged_data[,24650:24657])

test=apply(merged_data[,-c(1)], 2, scale)
head(test[,24650:24655])
test2=cbind.data.frame(merged_data[,1], test)
rownames(test2)=rownames(merged_data)
head(test2[,1:5])
merged_data_scaled=test2

#replace NA's with 0 (assume patients equal, average expression)
library(tidyr)
merged_data_scaled[is.na(merged_data_scaled)]=0
head(merged_data_scaled[,1:5])

colnames(merged_data_scaled)=c("severity_score", colnames(merged_data_scaled[,-c(1)]))
#calculate top correlated features now

merged_cor_scaled=cor(merged_data_scaled, use="pairwise.complete.obs")
summary(abs(merged_cor_scaled[,1]))
quantile(abs(merged_cor_scaled[,1]), 0.995, na.rm=T)

sig_cor_scaled=subset(merged_cor_scaled, abs(merged_cor_scaled[1,]) > 0.41)  #has the highest adjusted r-squared value
dim(sig_cor_scaled)
rownames(sig_cor_scaled)=c("severity_score", rownames(sig_cor_scaled[-c(1),]))
sig_cor_scaled[,1]
colnames(sig_cor_scaled)=c("severity_score", colnames(sig_cor_scaled[,-c(1)]))

scaled_features=subset(merged_data_scaled, select=rownames(sig_cor_scaled))
head(scaled_features)

fit_scaled=lm(severity_score ~., data=scaled_features)
summary(fit_scaled)

pred_r_squared(fit_scaled) #predicted r square = 0.4447


#look for genes that show clear patterns btwn high and low severity
head(sig_change_features, n=10)
head(change_high)
head(sort(change_high))

test=cbind.data.frame(merged_data[,1], merged_data$cd19_ENSG00000202538.1)
plot(test)
#all negatives gene exprs are above 2.5 severity, median 4.8
#updated prediction file with those test data w/ a negative gene expr



#look for other genes showing patterns in severity score
head(sort(change_low)[-c(1:20)])
test=cbind.data.frame(merged_data[,1], merged_data$cd8_ENSG00000202078.1)
plot(test)
#super high (>100) - above 4.5

test=cbind.data.frame(merged_data[,1], merged_data$cd8_ENSG00000279340.1)
plot(test)
#>50 range btwn 4 and 6.5

test=cbind.data.frame(merged_data[,1], merged_data$cd19_ENSG00000227470.1)
plot(test)
#when up, above 4.0

test=cbind.data.frame(merged_data[,1], merged_data$cd8_ENSG00000199824.1)
plot(test)
#highest one (huge outlier) had highest severity score

test=cbind.data.frame(merged_data[,1], merged_data$cd19_ENSG00000252614.1)
plot(test)
#above 100, above 4.0

test=cbind.data.frame(merged_data[,1], merged_data$nasal_ENSG00000163958.13)
plot(test)
#higher at least above 2.0

sig_cor_merged[1,]

test=cbind.data.frame(merged_data[,1], merged_data$cd4_ENSG00000136997.14)
plot(test)

test_fit=lm(severity_score~cd4_ENSG00000136997.14 + cd8_ENSG00000202078.1 + cd19_ENSG00000202538.1, data=merged_data_scaled)
test_fit$residuals

pred_r_squared(test_fit)
predict(test_fit, data.frame(test5))


library(umap)
umap=umap(test_features_nona, n_neighbors=7)
plot(umap$layout)
head(umap$layout)
head(merged_features)

umapplot = ggplot(data = data.frame(umap$layout), aes(x = X1, y = X2) ) + geom_point(aes(colour = scale(as.numeric(test_features_nona[,5]))), size = 5) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) + 
  
  #guides(colour=guide_legend(title=legendtitle)) +
  labs(colour = "Severity Score") + #scale_colour_distiller(palette = "Spectral")
  scale_colour_gradient(low="grey", high="red")

umapplot


head(test_features_nona[,1:5])
head(test_umap_features[,1:5])
test_umap_features=subset(test5, select=colnames(test_features_nona[,-c(1)]))
umap_test=umap(test_umap_features, n_neighbors=7)
plot(umap_test$layout)

umapplot = ggplot(data = data.frame(umap_test$layout), aes(x = X1, y = X2) ) + geom_point(aes(colour = scale(as.numeric(test_umap_features[,4]))), size = 5) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) + 
  
  #guides(colour=guide_legend(title=legendtitle)) +
  labs(colour = "Severity Score") + #scale_colour_distiller(palette = "Spectral")
  scale_colour_gradient(low="grey", high="red")

umapplot



#running lasso regression
library(caret)
library(glmnet)

ytrain=merged_data[,1]
xtrain=as.matrix(merged_data_scaled[,-c(1)])

lasso_cv=cv.glmnet(xtrain, ytrain, family="gaussian", alpha=1)
head(coef(lasso_cv), n=100)
plot(lasso_cv)
best_lambda=lasso_cv$lambda.min
cat(best_lambda)
lasso_mod=glmnet(xtrain, ytrain, family="gaussian", alpha=1, lambda=best_lambda)
coef(lasso_mod)

yhat=predict(lasso_mod, xtest)
head(yhat)

yhat

test_test=merged_data_scaled[10:58,]
test_train=merged_data_scaled[c(1:9,59:77),]

ytrain=test_train[,1]
xtrain=as.matrix(test_train[,-c(1)])

ytest=test_test[,1:2]
xtest=as.matrix(test_test[,-c(1)])

mse=mean((ytest-yhat)^2)
mse
yhat
ytest
result=cbind.data.frame(ytest[,1],yhat)
result


#lasso w/ features already deemed significant correlations
test_test=scaled_features[10:58,]
test_train=scaled_features[c(1:9,59:77),]

ytrain=test_train[,1]
xtrain=as.matrix(test_train[,-c(1)])

ytest=test_test[,1:2]
xtest=as.matrix(test_test[,-c(1)])

lasso_cv=cv.glmnet(xtrain, ytrain, family="gaussian", alpha=.5)
coef(lasso_cv)
plot(lasso_cv)
best_lambda=lasso_cv$lambda.min
cat(best_lambda)
lasso_mod=glmnet(xtrain, ytrain, family="gaussian", alpha=.5, lambda=best_lambda)
coef(lasso_mod)

yhat=predict(lasso_mod, xtest)
head(yhat)

yhat


mse=mean((ytest[,1]-yhat)^2)
mse
yhat
ytest
result=cbind.data.frame(ytest[,1],yhat)
result

plot=cbind.data.frame(scaled_features[,1], scaled_features$cd4_ENSG00000130638.15)
plot(plot)

plot1=cbind.data.frame(scaled_features[,1], scaled_features$cd8_ENSG00000099860.8)
plot(plot1)

plot2=cbind.data.frame(scaled_features[,1], scaled_features$cd8_ENSG00000164104.11)
plot(plot2)

plot3=cbind.data.frame(scaled_features[,1], scaled_features$nasal_ENSG00000005249.12)
plot(plot3)

plot4=cbind.data.frame(scaled_features[,1], scaled_features$nasal_ENSG00000148834.12)
plot(plot4)

plot5=cbind.data.frame(scaled_features[,1], scaled_features$nasal_ENSG00000174695.9)
plot(plot5)

lasso_test=lm(severity_score~cd4_ENSG00000130638.15+cd8_ENSG00000099860.8+cd8_ENSG00000164104.11+nasal_ENSG00000005249.12+nasal_ENSG00000148834.12+nasal_ENSG00000174695.9, data=scaled_features)
summary(lasso_test)
pred_r_squared(lasso_test)

save=predict(lasso_test, data.frame(merged_test_data_scaled))
save

predict=read.table("prediction.csv", sep=",")
head(predict)
predict=predict[-c(1),]
histogram(as.numeric(paste(predict[,1])), breaks=20)
as.numeric(paste(predict[,1]))
predict=predict[-c(1),]
median(as.numeric(paste(predict[,1])))
boxplot(as.numeric(paste(predict[,1])))
                  
quantile(as.numeric(paste(predict[,1])), 0.9)
