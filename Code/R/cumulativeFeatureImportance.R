library(ggplot2)

impindex = 3

rfmodel = readRDS("RDSFiles/rfmodel_comb.rds");
rankedt = (rfmodel$importance[order(-rfmodel$importance[,impindex]),]);

cname = rownames(rankedt)

nGram = 0;
nGDip = 0;
PSF   = 0;

# Uncomment if we take all features with positive importance score
# maxFeatures = min(which (rankedt[,impindex] < 0));
# cat("Count of features with positive importance:", maxFeatures,"\n");

# Uncomment when considering all features
maxFeatures = length(rankedt[,1]);

maxFeatures = 9936
for(i in 1:maxFeatures){
  if(gregexpr(pattern = "C_", cname[i])[[1]][1]>0){
    nGram = nGram + rankedt[i,impindex]
  } else if(gregexpr(pattern = "P_", cname[i])[[1]][1]>0){
    PSF = PSF + rankedt[i,impindex]
  } else {
    nGDip = nGDip + rankedt[i,impindex]
  }
}

cat("nGrams =", nGram,"\n")
cat("nGDip =" , nGDip,"\n")
cat("PSF ="   , PSF,"\n")

# Data for different feature types
data = data.frame(`Category` = character(), CumImp = numeric());
data = rbind(data, data.frame(`Category` = "PSN"    , CumImp = PSF));
data = rbind(data, data.frame(`Category` = "n-grams", CumImp = nGram));
data = rbind(data, data.frame(`Category` = "nGDip"  , CumImp = nGDip));

cumFeatureImp = ggplot(data,aes(x=`Category`, y = CumImp)) + 
  theme_bw(base_size = 40, base_family = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(aspect.ratio = 0.6) +
  geom_bar(aes(fill = Category),stat = "identity", position = "dodge", width = 0.5) +
  labs(x = "Feature Extraction Technique", y = "Cumulative Importance") +
  ylim(0,0.15);
  
postscript(file = paste0("cumFeatureImp", maxFeatures, ".eps"), paper = "letter");
cumFeatureImp;
dev.off();