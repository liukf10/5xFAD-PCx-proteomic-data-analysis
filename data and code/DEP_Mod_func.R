FCS <- function(stat ,  moduleinf, DEP = NULL, ordername = c("fc","p"),filter = "p",
                fc = 1.2, fcname = "fc", p = 0.05, pname = c("p","p.adj"),...){
  if(!missing(stat)) {
    if(length(ordername) > 1) ordername = match.arg(ordername, c("fc","p"));
    if(length(pname) > 1) pname = match.arg(pname,c("p","p.adj"));
    if(!all(c(ordername,pname,fcname) %in% colnames(stat))) stop("wrong ordrename or fcname")
    stat <- stat[!is.na(stat[,pname]),];
    DEP <- stat[stat[,fcname] >= fc | stat[,fcname] <= 1/fc,];
    DEP <- DEP[DEP[,pname] <= p ,];
    if(ordername == "fc") DEP <- DEP[order(abs(log2(DEP[,fcname])), decreasing = T),];
    if(ordername == "p") DEP <- DEP[order(DEP[,"p"]),];
  }
  if(!"ori.ID" %in% colnames(DEP)) stop("No DEP ori.ID.")
  FCSen <- Module_Enrich(moduleinf, DEP$ori.ID, coln="ori.ID");
  Z <- -log10(FCSen$p.adj);
  Zmean <- apply(Z,2,mean);
  num1 <- apply(FCSen$Counts,2,function(x) which(x!=0)[1]);
  num1[is.na(num1)] <- 0;
  Zmean <- Zmean*nrow(Z)/(nrow(Z)-num1);
  pmean <- 10^-Zmean;
  FCSenrich <- FCSenrichplot(FCSen, filter = filter,...);
  list(DEPinf = DEP, FCSenrichinf= FCSen,FCSplotinf = FCSenrich,pmean,Zmean)
}



#DEPmod enrich heatmap
enrichmod <- function(Module, mod, ID, IDtype = "all", coln = "ori.ID", 
                      datainf = NULL, FCS = TRUE, 
                      padj = "fdr", filter = "p.adj", removeMod = NULL,
                      filename = NULL, filetype = "pdf") {
  IDtype = match.arg(IDtype, c("all","up","down","none"));
  if (missing(mod))  IDtype = "none"
  if (IDtype == "all") ID = mod$DEPinf$ori.ID;
  if (IDtype == "up") ID = mod$DEPinf$ori.ID[mod$DEPinf$fc>1];
  if (IDtype == "down") ID = mod$DEPinf$ori.ID[mod$DEPinf$fc<1];
  ORAenrich <- Module_Enrich(Module, ID, coln = coln, 
                             enrichtype = "ORA", datainf = datainf,
                             p.adj.method = padj)
  rowname <- ORAenrich$module.name;
  ORAenrich <- data.frame(Counts = ORAenrich$Counts,
                          module.size = ORAenrich$module.size,
                          precent = ORAenrich$precent,
                          p = ORAenrich$p, p.adj = ORAenrich$p.adj,
                          Z.score = ORAenrich$Z.score,
                          stringsAsFactors = FALSE)
  rownames(ORAenrich) <- rowname;
  rm(rowname)
  FCSen <- Module_Enrich(Module, ID, coln=coln);
  Z <- -log10(FCSen$p.adj);
  Zmean <- apply(Z,2,mean);
  Zmean2 <- apply(FCSen$Z.score,2,mean);
  num1 <- apply(FCSen$Counts,2,function(x) which(x!=0)[1]);
  num1[is.na(num1)] <- 0;
  Zmean <- Zmean*nrow(Z)/(nrow(Z)-num1); Zmean2 <- Zmean2*nrow(Z)/(nrow(Z)-num1);
  pmean <- 10^-Zmean; pmean2 <- 10^-Zmean2;
  ORAenrich$p.adj <- pmean; ORAenrich$Z.score <- Zmean; ORAenrich$p <- pmean2;
  if (FCS) {
    FCSenrich <- Module_Enrich(Module, ID, coln = coln, datainf = datainf,
                               p.adj.method = padj)
    if (filter == "none" && !is.null(removeMod)) {
      a <- removeMod < 10 & removeMod > 0;
      b <- !a;
      removeMod[a] <- paste0("M0",removeMod[a]);
      removeMod[b] <- paste0("M",removeMod[b]);
      rm(a,b);
      remove <- FCSenrich$module.name %in% removeMod;
      FCSenrich$Counts <- FCSenrich$Counts[,!remove]
      FCSenrich$module.size <- FCSenrich$module.size[!remove]
      FCSenrich$module.name <- FCSenrich$module.name[!remove]
      FCSenrich$precent <- FCSenrich$precent[,!remove]
      FCSenrich$p <- FCSenrich$p[,!remove]
      FCSenrich$p.adj <- FCSenrich$p.adj[,!remove]
      FCSenrich$Z.score <- FCSenrich$Z.score[,!remove]
    }
    FCSenrich <- FCSenrichplot(FCSenrich, filter = filter,
                               filename = filename, filetype = filetype)
    mod <- levels(as.factor(FCSenrich$module))
    enrich <- ORAenrich[rownames(ORAenrich) %in% mod, ]
  } else ORAenrich
}

if(1!=1){
  DEP_Mod_HeatMap <- function(DEP_Mod, filter = c("p","p.adj"), 
                              cutoff = 0.05, filename = NULL, ...) {
    filter <- match.arg(filter, c("p","p.adj"));
    if (!is.list(DEP_Mod)) stop("DEP_Mod is not the list.")
    x <- matrix(data = NA, nrow = length(DEP_Mod[[1]][[1]]), ncol = length(DEP_Mod))
    x <- as.data.frame(x);  rownames(x) <- rownames(DEP_Mod[[1]])
    if (is.null(names(DEP_Mod))) names(DEP_Mod) <- 1:length(DEP_Mod);
    colnames(x) <- names(DEP_Mod)
    x1 = x2 = x3 = x4 = x5 = x
    for (i in 1:length(DEP_Mod) ) {
      if (all(c("precent","Counts","module.size","module.size","p","p.adj") %in% 
              names(DEP_Mod[[i]]) ) ) {
        if(length(DEP_Mod[[i]][["precent"]]) != nrow(x1))
          stop ( paste(names(DEP_Mod)[i], "is not the same module number with others.") );
        x1[ , i ] = DEP_Mod[[i]][["precent"]];
        x2[ , i ] = DEP_Mod[[i]][["Counts"]];
        x3[ , i ] = DEP_Mod[[i]][["module.size"]];
        x4[ , i ] = DEP_Mod[[i]][["p"]];
        x5[ , i ] = DEP_Mod[[i]][["p.adj"]];
      } else stop( paste(names(DEP_Mod)[i], "is not the correct data."))
    }
    if (filter == "p")
      connect <- list(precent = x1, Counts =x2,module.size = x3, p = x4)
    if (filter == "p.adj")
      connect <- list(precent = x1, Counts =x2,module.size = x3, p = x5)
    p <- connect$p
    p <- signif(p, digits = 2)
    p <- format(p, scientific = TRUE)
    p[connect$p > cutoff] <- ""
    Counts <- connect$Counts
    Counts[connect$p > cutoff] <- ""
    textMatrix = paste(as.matrix(Counts),"\n(",as.matrix(p),")",sep="")
    dim(textMatrix) = dim(Counts)
    textMatrix[textMatrix=="\n()"]<-""
    precent <- connect$precent
    ratio <- colSums(connect$Counts) / colSums(connect$module.size)
    enrichFold <- t(t(precent)/ratio/100)
    if (!is.null(filename)) pdf(paste0("plot/",filename,".pdf"))
    if(1!=1)
      WGCNA::labeledHeatmap(Matrix = enrichFold, xLabels = colnames(p), 
                            yLabels = rownames(p), 
                            cex.lab = 1, colorLabels = TRUE, 
                            colors = WGCNA::blueWhiteRed(100)[51:100], 
                            textMatrix = textMatrix, 
                            setStdMargins = FALSE, 
                            cex.text = 1, ...)
    if(1==1)
      WGCNA::labeledHeatmap(Matrix = t(enrichFold), xLabels = rownames(p), 
                            yLabels = colnames(p), 
                            cex.lab = 1, colorLabels = TRUE, 
                            colors = WGCNA::blueWhiteRed(100)[51:100], 
                            textMatrix = t(textMatrix), 
                            setStdMargins = FALSE, 
                            cex.text = 1, ...)
    if (!is.null(filename)) dev.off()
    list(enrichFold = enrichFold,
         textMatrix = textMatrix)
  }
}
