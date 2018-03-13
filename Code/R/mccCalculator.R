# data = read.csv("temp.csv")

Sn  = c(.9250, .8279, .6237)
Sp  = c(.6560, .5913, .7741	)
Acc = c(.7900, .7096, .6989	)

# Sn  = data$Sensitivity
# Sp  = data$Specificity
# Acc = data$Accuracy

P   = 93;
N   = 93;
# P = 73;
# N = 1390;

mcc =c();
prec = c();


for (i in 1:length(Sn)) {
  TP = Sn[i] * P;
  TN = Sp[i] * N;
  
  FP = N - TN;
  FN = P - TP;
  
  curMcc = (TP * TN - FP * FN)/sqrt((TP + FP)*(TP + FN)*(TN + FN)*(TN + FP))
  curPrec = TP / (TP + FP)
  
  cat("DiffCheck:", TP + TN - Acc[i] * (P + N), "mcc:", curMcc, "prec:", curPrec, "\n")
  
  mcc  = c(mcc, curMcc);
  prec = c(prec, curPrec);
}
  
