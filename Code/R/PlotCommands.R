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
