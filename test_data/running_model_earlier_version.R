
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
