DAPextract <- function(data, moduleinf, moduleNum, markedID ,
                       IDname = NULL, coln = "new.ID",OnlyPlotLast = TRUE,...) {
  Mod <- getmoduleHub(data, moduleinf, moduleNum, coln = coln, adjustp = FALSE,...);
  IDsets = list(gene.set = list(ID = markedID), color.code = "red");
  names(IDsets$gene.set) <- deparse(substitute(markedID));
  if (!is.null(IDname)) names(IDsets$gene.set) <- IDname;
  DEPnet <- try(DEP_Mod_net_plot(Mod, IDsets = IDsets, data = data, module = moduleinf,
                                 filename = NULL, reconstructNet = TRUE,OnlyPlotLast = OnlyPlotLast,...), silent = T);
  if ( class(DEPnet) == "try-error") 
    DEPnet <- DEP_Mod_net_plot(Mod, IDsets = IDsets, data = data, module = moduleinf, 
                               filename = NULL, reconstructNet = F,OnlyPlotLast = OnlyPlotLast,...);
  DAPID <- DEPnet$netgene;
}
DAPbyMod <- function(data, moduleinf, moduleNum, markedID,
                     IDname = NULL, coln = "new.ID",...) {
  DAP = NULL;
  for( i in 1:length(moduleNum)) {
    DAPID <- DAPextract(data, moduleinf, moduleNum[i], markedID, IDname = IDname, coln = coln,...);
    names(DAPID) <- rep(moduleNum[i],length(DAPID));
    DAP <- c(DAP, DAPID);
  }
  DAP
}


DAPextract2 <- function(data, moduleinf, moduleNum, markedID ,
                        IDname = NULL, coln = "new.ID",OnlyPlotLast = TRUE,...) {
  Mod <- getmoduleHub(data, moduleinf, moduleNum, coln = coln, adjustp = FALSE,...);
  IDsets = list(gene.set = list(ID = markedID), color.code = "red");
  names(IDsets$gene.set) <- deparse(substitute(markedID));
  if (!is.null(IDname)) names(IDsets$gene.set) <- IDname;
  DEPnet <- try(DEP_Mod_net_plot(Mod, IDsets = IDsets, data = data, module = moduleinf,
                                 filename = NULL, reconstructNet = TRUE,OnlyPlotLast = OnlyPlotLast,...), silent = T);
  if ( class(DEPnet) == "try-error") 
    DEPnet <- DEP_Mod_net_plot(Mod, IDsets = IDsets, data = data, module = moduleinf, 
                               filename = NULL, reconstructNet = F,OnlyPlotLast = OnlyPlotLast,...);
  DEPnet;
}
