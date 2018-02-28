#########################################################
# For writing to R window comment out the following line
# in each graph
# theme_bw(base_size = 36, base_family = "") +
#########################################################

library("ggplot2")
library("XLConnect")
library("reshape2")

###### Accuracy/MCC etc. vs. choice of nFeatures  ############

# Use the appropriate data file here:
xlsFile  = "PerfSearch_comb.xlsx"
xlsSheet = "TenfoldCV_SVM_Micro"

workBook = loadWorkbook(xlsFile)
data = readWorksheet(workBook, xlsSheet);

data = data[, c("nF", "AUCROC", "Accuracy", "Sensitivity", "Specificity", "MCC")]
colnames(data)[2] = "auROC"

df <- melt(data,  id.vars = "nF", variable.name = 'Metric');

nFeatureTuning = ggplot(df,aes(x=nF,y=value)) +
  theme_bw(base_size = 36, base_family = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top") +
  geom_line(aes(colour=Metric),size =3) +
  labs(x = "Number of Features", y = "Performance Score x 100");
  
postscript(file = paste0(xlsSheet, ".eps"), paper = "letter");
nFeatureTuning;
dev.off();

###### feature Space Perf Comparisons ###################

# Data for different feature types within top 100 features
epsFile = "featureSpacePerfCmp.eps";
data = NULL;
data = rbind(data, c(0.55, 0.09, 1.00, 0.21));
data = rbind(data, c(0.65, 0.82, 0.48, 0.32));
data = rbind(data, c(0.79, 0.76, 0.81, 0.57));
data = rbind(data, c(0.87, 0.83, 0.91, 0.74));
rownames(data) = c("PSF", "n-grams", "nGDip", "COMB");
colnames(data) = c("Accuracy", "Sensitivity", "Specificity", "MCC");

# Data for different feature types each using 100 features
epsFile = "featureSpacePerfCmp2.eps";
data = NULL;
data = rbind(data, c(0.62, 0.28, 0.95, 0.32));
data = rbind(data, c(0.74, 0.80, 0.68, 0.49));
data = rbind(data, c(0.84, 0.84, 0.83, 0.68));
data = rbind(data, c(0.87, 0.83, 0.91, 0.74));
rownames(data) = c("PSF", "n-grams", "nGDip", "COMB");
colnames(data) = c("Accuracy", "Sensitivity", "Specificity", "MCC");

# Data for subset combinations of top 100 features
epsFile = "featureSpacePerfCmp3.eps";
data = NULL;
data = rbind(data, c(0.84, 0.81, 0.87, 0.68));
data = rbind(data, c(0.73, 0.74, 0.71, 0.46));
data = rbind(data, c(0.80, 0.75, 0.85, 0.60));
data = rbind(data, c(0.87, 0.83, 0.91, 0.74));
  # COMB1 = nGrams(38) + nGDip(62)
  # COMB2 = nGrams(67) + PSF(33)
  # COMB3 = nGDip(75) + PSF(25)
rownames(data) = c("COMB1", "COMB2", "COMB3", "COMB");
colnames(data) = c("Accuracy", "Sensitivity", "Specificity", "MCC");

# Draw the graph
df <- melt(data);
featureSpacePerfCmp = ggplot(df,aes(x=Var1, y = value)) +
  theme_bw(base_size = 36, base_family = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top") +
  theme(aspect.ratio = 0.6) +
  geom_bar(aes(fill = Var2),stat = "identity", position = "dodge") +
  labs(x = "Feature Extraction Technique", y = "Performance score") +
  ylim(0,1)

postscript(file = epsFile, paper = "letter");
featureSpacePerfCmp;
dev.off();

# Feature weights bar plot

svmModel = readRDS("svmDNAProt.rds");
w = t(svmModel$coefs) %*% svmModel$SV

# sort by absolute value
w = w[,order(-abs(w))]

df = melt(w[1:25])

featureWeightsBar = ggplot(df,aes(x=reorder(rownames(df), -value), y = value)) +
  theme_bw(base_size = 18, base_family = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 0.6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_bar(aes(fill=(value > 0)),stat = "identity", position = "dodge") +
  labs(x = "Features", y = "Discriminant weight")
 
 postscript("featureWeightsBar.eps", paper = "letter");
 featureWeightsBar;
 dev.off();

