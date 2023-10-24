library(corrplot)
library("Hmisc")

df_f <-read.csv("Data/Derived/master_f.csv")
df_uf <- read.csv("Data/Derived/master_uf.csv")

#taking onlu the vlaues 
df_f$D_change <- df_f$D_00 - df_f$D_51
cor_master <- df_uf[,c(5:9, 16,19,17,18,32, 20:26)]
cor_master <- na.omit(cor_master)

cor_master.cor = cor(cor_master)

cor_master.cor = cor(cor_master, method = c("spearman"))

#look at p values of correlations 
cor_master.rcor = rcorr(as.matrix(cor_master))
cor_master.coeff = cor_master.rcorr$r
cor_master.p = cor_master.rcorr$p

corrplot(cor_master.cor)
