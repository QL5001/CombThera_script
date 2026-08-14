#Principal component analysis (PCA) of RNA-Seq reads of Enterococcus faecalis OG1RF strain for determining altered transcriptomic changes under the treatment of Streptomyces-dervied cyclic short-chain antimicrobial peptide siamycin.#
#Begin#
library(knitr)
library(mixOmics)
library(dplyr)
library(stats)

RS_Sia1 = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)    #open file "group_S_Differential_Expression_PCA" consisting the differentially-regulated E. faecalis OG1RF transcriptomes under siamycin treatment

X = RS_Sia1[, 1:2414]  
X_mean = apply(X, 2, mean)
X_center = scale(X, center = X_mean, scale = FALSE)
Y = RS_Sia1$Sample

pca.RS = pca(X_center, ncomp = 6, center = FALSE, scale = FALSE)
plot(pca.RS, ylim = c(0, 1.0), main = "PCA on RS")

ev = pca.RS$prop_expl_var   #explained_variance
print(ev)
PCA_ev = data.frame(ev)
write.csv(PCA_ev, "PCA_ev_241011.csv")  #export explained variance for each principal component into a csv file

sco = pca.RS$variates   #prinicipal components values
print(sco)
sco_PC = data.frame(SampleNo = RS_Sia1$Sample, sco)
write.csv(sco_PC, "PCA_PC_241011.csv")   #export the principal component values into a csv file

#Exported csv files associated with principal components and explained variances will be used for further analyses and creating figures in Excel, Prism, and R programming
#End#
