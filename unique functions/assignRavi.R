assignRavi <- function (object, RG.features,  RH.features, RI.features, RN.features ,set.ident = FALSE, ctrl_genes = 100,
                            ...)
{
  name <- "Cell.Cycle"
  features <- list(  RG.Score = RG.features, RH.Score =  RH.features, RI.Score = RI.features, RN.Score = RN.features)
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
                                                                 first = "Radial Glia", second = "Reactive Hypoxia", third = "Reactive Immune", fourth = "Regional NPC", null = "Undecided") {
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
  colnames(x = cc.scores) <- c("rownames", "RG.Score", "RH.Score", "RI.Score","RN.Score",
                               "Ravi")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("RG.Score", "RH.Score", "RI.Score" ,"RN.Score","Ravi")]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Ravi"
  }
  return(object)
}

