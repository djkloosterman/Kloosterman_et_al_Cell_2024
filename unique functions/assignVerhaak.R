assignVerhaak <- function (object, MES.features,  PN.features, CL.features ,set.ident = FALSE, ctrl_genes = 100,
                           ...)
{
  name <- "Cell.Subtype"
  features <- list(  MES.Score = MES.features, PN.Score =  PN.features, CL.Score = CL.features)
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
                                                                 first = "MES", second = "PN", third = "CL",  null = "Undecided") {
    if (all(scores < -0)) {
      return(null)
    }
    else {
      return(c(first, second, third)[which(x = scores == max(scores))])
    }
  }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments),
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "MES.Score", "PN.Score", "CL.Score",
                               "Verhaak")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("MES.Score", "PN.Score", "CL.Score","Verhaak")]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Verhaak"
  }
  return(object)
}
