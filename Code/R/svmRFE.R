library("ROCR")

  svmRFE <-
    function(features, labels, jump) {
      nF = ncol(features);
      
      indSurvivor = seq(1:nF);
      rankedFeatures = vector(length=nF);
      indRankedFeature = nF;

      while(length(indSurvivor) > 0) {
        cat(as.character(Sys.time()),">> Remaining feature count:", length(indSurvivor), "\n");
        
        svmModel = svm(features[, indSurvivor], labels, kernel="linear");
        
        svmpred = predict(svmModel, features[, indSurvivor])
        svmprediction = prediction(as.numeric(svmpred), as.numeric(labels))
        
        nSV  = svmModel$tot.nSV;
        acc  = round(unlist(ROCR::performance(svmprediction,"acc")@y.values)[2], 2);
        sens = round(unlist(ROCR::performance(svmprediction,"sens")@y.values)[2], 2);
        spec = round(unlist(ROCR::performance(svmprediction,"spec")@y.values)[2], 2);
        prec = round(unlist(ROCR::performance(svmprediction,"prec")@y.values)[2], 2);
        mcc  = round(unlist(ROCR::performance(svmprediction,"mat")@y.values)[2], 2);
        cat(as.character(Sys.time()),">> nSV:", nSV, ", Acc: ", acc, ", Sens: ", sens,
            ", Spec: ", spec, ", Prec: ", prec, ", MCC: ", mcc, "\n");
        
        #compute the weight vector
        w = t(svmModel$coefs) %*% svmModel$SV
        
        #compute ranking criteria
        rankingCriteria = w * w
        
        #rank the features
        ranking = order(rankingCriteria);
        
        for (i in 1:min(jump,length(ranking))) {
          #update feature ranked list
          rankedFeatures[indRankedFeature] = indSurvivor[ranking[i]]
          indRankedFeature = indRankedFeature - 1
        }
        
        #eliminate the feature with smallest ranking criterion
        indSurvivor = indSurvivor[-ranking[1:jump]]
        
        write.csv(rankedFeatures, "rankedFeatures.csv");
        write.csv(indSurvivor, "indSurvivor.csv");
      }
      
      # We have so far maintained the indices of the features
      # Now convert to the actual feature names
      rankedFeatures = colnames(features)[rankedFeatures];
      
      return (rankedFeatures)
    }