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

