assignLocation <- function (object, ct.features,  pan.features, le.features, mvp.features ,set.ident = FALSE, ctrl_genes = 100,
                           ...)
{
  name <- "Cell.Cycle"
  features <- list(  PAN.Score = pan.features, LE.Score =  le.features, MVP.Score = mvp.features, CT.Score = ct.features)
  object.cc <- AddModuleScore(object = object,
                              features = features,
                              name = name,
                              ctrl = min(c(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1))), 
                                         ctrl_genes), 
                              ...)
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]),
                     value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  rm(object.cc)
  CheckGC()
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores,
                                                                 first = "PAN", second = "LE", third = "MVP", fourth = "CT", null = "Undecided") {
    if (all(scores < -0)) {
      return(null)
    }
    else {
      return(c(first, second, third, fourth)[which(x = scores == max(scores))])
    }
  }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments),
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "PAN.Score", "LE.Score", "MVP.Score","CT.Score",
                               "Location")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("PAN.Score", "LE.Score", "MVP.Score" ,"CT.Score","Location")]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Location"
  }
  return(object)
}
