filteredFeaturization <-
  function(sequences, featureNames) {
    features = data.frame(1:length(sequences))
    # a dummy column got created. Let us name it. We will
    # delete this column at the end
    colnames(features)[length(features)] = "Serial"
    
    gapLengths = NULL;
    nMerLen = 0;
    
    featureNameMap = new.env();
    for (key in featureNames) {
      # create the column
      features[key] = integer(nrow(features));
      # save the feature name in a map
      assign(key, TRUE, featureNameMap);
  
      splits = strsplit(toString(key), "_")[[1]];
      type = splits[1];
      if (splits[1] == 'G') {
        gap = as.numeric(splits[2]);
        if (!(gap %in% gapLengths)) {
          gapLengths = c(gapLengths, gap);
        }
      } else {
        nMerLen = max(nMerLen, nchar(splits[3]));
      }
    }
    if (!is.null(gapLengths)) {
      gapLengths = gapLengths[order(gapLengths)];
    }
    
    for (i in 1:nrow(features)) {
      strSeq = strsplit(toString(sequences[i]), "")[[1]];
      for (j in 1:length(strSeq)) {
        token = "";
        for (k in 1:nMerLen) {
          if (j+k-1 > length(strSeq)) {
            break;
          }

          token = paste(token, strSeq[j+k-1], sep = "");
          
          # update the n-mer feature count
          countToken = paste("C", 0, token, sep = "_")
          if (exists(countToken, envir = featureNameMap)) {
            features[i,countToken] = features[i,countToken] + 1;
          }
          
          # update the posorder feature count
          posToken = paste("P", j, token, sep = "_");
          if (j <= 10 && exists(posToken, envir = featureNameMap)) {
            features[i,posToken] = features[i,posToken] + 1;
          }
        }
      }
    }

    # gapped dipeptide feature counts
    if (!is.null(gapLengths)) {
      for (i in 1:nrow(features)) {
        strSeq = strsplit(toString(sequences[i]), "")[[1]];
        for (j in 1:length(strSeq)) {
  
          for (k in gapLengths) {
            if (j+1+k > length(strSeq)) {
              break;
            }
  
            token = paste(strSeq[j], strSeq[j+1+k], sep = "");
            token = paste("G", k, token, sep = "_");
            if (exists(token, envir = featureNameMap)) {
              features[i,token] = features[i,token] + 1;
            }
          }
        }
      }
    }
    features$Serial = NULL
    
    return(features)
  }