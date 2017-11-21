#' Illustration of Featurization
#'
#' This function takes dataset and extracts features.
#'
#' @param sequences provided as dataframe
#' @param labels class label of each data. Provided as a column vector
#' @param alphabet is the list of aminoacids
#'        c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y")
#' @param nGramOrder Highest value of n in n-grams feature extraction technique
#' @param nGDipOrder Highest value of n in n-Gapped-Dipeptide (nGDip) feature extraction technique
#' @param psfOrder Highest value of n in n-grams in Position Specific Feature (PSF) feature extraction technique
#' @return a featurized dataframe
#' @export
featurization <-
  function(sequences, labels, alphabet, nGramOrder, nGDipOrder, psfOrder) {
    features = data.frame(1:length(sequences))
    # a dummy column got created. Let us name it. We will
    # delete this column at the end
    colnames(features)[length(features)] = "Serial"
    
    alphaMap = new.env();
    for (key in alphabet) {
      assign(key, TRUE, alphaMap);
    }
    
    nGramCount = 0;
    psfCount = 0;
    nGDipCount = 0;
    
    if (nGramOrder > 0 || psfOrder > 0) {
      for (i in 1:nrow(features)) {
        strSeq = strsplit(toString(sequences[i]), "")[[1]];
        for (j in 1:length(strSeq)) {
          token = "";
          for (k in 1:max(nGramOrder, psfOrder)) {
            if (j+k-1 > length(strSeq)) {
              break;
            }
            if (!(exists(strSeq[j+k-1], envir = alphaMap))) {
              break;
            }

            token = paste(token, strSeq[j+k-1],sep = "");
            
            # update the nGramOrder feature count
            if (nchar(token) <= nGramOrder) {
              countToken = paste("C", 0, token, sep = "_")
              if (!(countToken %in% colnames(features))) {
                # create the column on demand
                features[countToken] = integer(nrow(features));
                nGramCount = nGramCount + 1;
              }
              features[i,countToken] = features[i,countToken] + 1/(length(strSeq) - k + 1);
            }
            
            # update the psfOrder feature count
            if (j <= 10 && nchar(token) <= psfOrder) {
              posToken = paste("P", j, token, sep = "_");
              if (!(posToken %in% colnames(features))) {
                # create the column on demand
                features[posToken] = integer(nrow(features));
                psfCount = psfCount + 1;
              }
              features[i,posToken] = features[i,posToken] + 1;
            }
          }
        }
      }
    }
    cat(as.character(Sys.time()),">> n-grams based features:", nGramCount, "\n");
    cat(as.character(Sys.time()),">> Position Specific Features (PSF):", psfCount, "\n");
    
    if (nGDipOrder > 0) {
  
      for (i in 1:nrow(features)) {
        strSeq = strsplit(toString(sequences[i]), "")[[1]];
        for (j in 1:length(strSeq)) {
          if (!(exists(strSeq[j], envir = alphaMap))) {
            next;
          }

          for (k in 1:nGDipOrder) {
            if (j+1+k > length(strSeq)) {
              break;
            }
            if (!(exists(strSeq[j+1+k], envir = alphaMap))) {
              next;
            }
            token = paste(strSeq[j], strSeq[j+1+k], sep = "");
            token = paste("G", k, token, sep = "_");
            if (!(token %in% colnames(features))) {
              # create the column on demand
              features[token] = integer(nrow(features));
              nGDipCount = nGDipCount + 1;
            }
            
            features[i,token] = features[i,token] + 1/(length(strSeq) - k -1);
          }
        }
      }
    }
    cat(as.character(Sys.time()),">> n-Gapped-Dipeptide (nGDip) based features:", nGDipCount, "\n");
    
    features$protection = as.factor(labels);
    names(features) = make.names(names(features));
    features$Serial = NULL;
    cat(as.character(Sys.time()),">> Total features: ", length(features[1,]) - 1, "\n");
    
    return(features)
  }