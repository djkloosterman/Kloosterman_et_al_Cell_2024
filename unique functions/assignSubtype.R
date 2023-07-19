assignSubtype <- function (object, MES1.features,  MES2.features, AC.features, OPC.features, NPC1.features, NPC2.features ,set.ident = FALSE, ctrl_genes = 100,
                           ...)
{
  name <- "Cell.Subtype"
  features <- list(  MES1.Score = MES1.features, MES2.Score =  MES2.features, AC.Score = AC.features, OPC.Score = OPC.features, NPC1.Score = NPC1.features, NPC2.Score = NPC2.features)
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
                                                                 first = "MES1", second = "MES2", third = "AC", fourth = "OPC", fifth = "NPC1", sixth = "NPC2", null = "Undecided") {
    if (all(scores < -0)) {
      return(null)
    }
    else {
      return(c(first, second, third, fourth, fifth, sixth)[which(x = scores == max(scores))])
    }
  }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments),
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "MES1.Score", "MES2.Score", "AC.Score","OPC.Score", "NPC1.Score", "NPC2.Score",
                               "Cell.Subtype")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("MES1.Score", "MES2.Score", "AC.Score","OPC.Score", "NPC1.Score", "NPC2.Score","Cell.Subtype")]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Cell.Subtype"
  }
  return(object)
}
