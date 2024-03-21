message("rrho_data func need statinf matrix which included idname, pname and logFCname column.")
rrho_data <- function(dat, idname = "ID", pname = "p.value",
                      logFCname = "log2FC") {
  if(!idname %in% colnames(dat)) stop("wrong idname.");
  if(!pname %in% colnames(dat)) stop("wrong pname.");
  if(!logFCname %in% colnames(dat)) stop("wrong logFCname.");
  dat2 <- data.frame(ID = dat[,idname], Z = -log10(dat[,pname]));
  dat2$Z[dat[,logFCname] < 0] <- -dat2$Z[dat[,logFCname] < 0];
  dat2
}

rrho_web <- function(rrho_data1, rrho_data2) {
  if(!"ID" %in% colnames(rrho_data1)) stop("rrho_data1 should have column: ID.");
  if(!"Z" %in% colnames(rrho_data1)) stop("rrho_data1 should have column: Z.");
  if(!"ID" %in% colnames(rrho_data2)) stop("rrho_data2 should have column: ID.");
  if(!"Z" %in% colnames(rrho_data2)) stop("rrho_data2 should have column: Z.");
  pos = match(rrho_data1$ID,rrho_data2$ID);
  n1 <- nrow(rrho_data1); n2 <- nrow(rrho_data2);
  message(paste0("The original gene number is ",n1,", ",n2))
  rrho_data1 <- rrho_data1[which(!is.na(pos)),];
  rrho_data2 <- rrho_data2[pos[!is.na(pos)],];
  n3 <- nrow(rrho_data1);
  message(paste0("The overlap gene number is ",n3));
  if (n3 != 0) {
    rrho_data1$rank1 <- rank(rrho_data1$Z);
    rrho_data2$rank2 <- rank(rrho_data2$Z);
    rrho <- data.frame(ID1 = rrho_data1$ID, ID2 = rrho_data2$ID,
                      rank1 = rrho_data1$rank1, rank2 = rrho_data2$rank2,
                      metric1 = rrho_data1$Z, metric2 = rrho_data2$Z)
    } else stop ("No overlaped gene in two data.")
  rrho}


rrho_plot <- function(rrho_data, labels, 
                      color = NA, color_range = NA,
                      file = NA) {
  if(is.na(color)) {
    jet.colors <- colorRampPalette(c("#00007F", 
                                     "blue", "#007FFF", "cyan", 
                                     "#7FFF7F", "yellow", "#FF7F00", 
                                     "red", "#7F0000"))
    color = jet.colors(100)
  }
  if(missing(labels)) {stop("no labels.")}
  x <- data.frame(t(rrho_data$hypermat));
  maxn <- max(nchar(colnames(x)));
  while(any(nchar(colnames(x)) < maxn)) {
    colnames(x)[nchar(colnames(x))<maxn] <- gsub("X","X0",colnames(x)[nchar(colnames(x))<maxn])
  }
  maxn <- max(nchar(rownames(x)));
  while(any(nchar(rownames(x)) < maxn)) {
    rownames(x)[nchar(rownames(x))<maxn] <- paste0(0,rownames(x)[nchar(rownames(x))<maxn])
  }
  if(is.na(color_range[1])) color_range = c(min(x),max(x));
  p <- ggheatmap(x,color = color,
                 text_show_rows = "a", text_show_cols = "b",
                 text_position_rows = "left",
                 legendName = "-logP") +
    scale_fill_gradientn(colours = color,limit = color_range)+
    xlab(labels[1])+ylab(labels[2])+
    theme(axis.title.x = element_text(), axis.title.y = element_text())
  if(!is.na(file)) 
    jpeg(file = paste0("RRHO/",file,"_MAP_",labels[1],"_vs_",labels[2],".jpg"),
        width = 8, height = 8, units = "in", quality = 100, res = 300)
  print(p)
  if(!is.na(file)) graphics.off()
  p
}
