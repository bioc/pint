

join.top.regions <- function (model, feature.info, quantile.th = 0.95, augment = FALSE) {

  # Pick all top regions (based on the given threshold) and join overlapping windows
  # to get easily interpretable summaries

  # model: pint model object
  # feature.info: for instance ge$info; location info for genes
  # ntop: number top windows to check

  feature.info <- order.feature.info(feature.info)
  df <- as.data.frame(model)
  th <- quantile(df[["dependencyScore"]], quantile.th)
  #topg <- subset(df, dependencyScore > th)$geneName
  topg <- df$geneName[df$dependencyScore > th]

  # Get window around each top gene
  mygenes <- lapply(topg, function(gn) { rownames(findModel(model, gn)@W$X) })

  # Add other known genes in each region 
  if (augment && !is.null(feature.info)) {
    mygenes.aug <- lapply(mygenes, function (reg) { augment.region(reg, feature.info) })
  }

  # List all genes residing on the top regions
  topw <- unique(unname(unlist(mygenes.aug)))

  # Match to feature table, which is ordered(!) by chromosomal locations
  ord <- match(topw, rownames(feature.info))
  o <- order(ord) 
  ord <- ord[o]

  # Mark break points
  do <- c(1, diff(ord) == 1) 

  regs <- list()
  i <- 1
  k <- 0
  reg <- c()
  while (i <= length(do)) {
  
    if (do[[i]] == 1) {
      reg <- c(reg, ord[[i]])
    } else {
      k <- k + 1
      regs[[k]] <- reg
      reg <- c(ord[[i]])
    }
    i <- i + 1
  }

  regs <- lapply(regs, function (x) {rownames(feature.info)[x]})
  names(regs) <- as.character(1:length(regs))

  regs

}


summarize.region.parameters <- function (region.genes, model, X, Y) {

  # Take average of the Zs and Ws over the overlapping models
  # Useful for interpretation.
  
  zs <- array(NA, dim = c(ncol(X$data), length(region.genes)))
  wxs <- wys <- array(NA, dim = c(length(region.genes), length(region.genes)))
  colnames(zs) <- colnames(wxs) <- colnames(wys) <- rownames(wxs) <- rownames(wys) <- region.genes
  rownames(zs) <- colnames(X$data)

  chr <- na.omit(as.numeric(as.character(unique(X$info[region.genes,"chr"]))))
  chr.models <- model@chromosomeModels[[chr]]@models
  # FIXME: provide function that automatically fetches model names!
  model.names <- sapply(1:length(chr.models), function (k) {chr.models[[k]]@geneName})
  # Pick genes from this region that also have models
  gs <- intersect(model.names, region.genes)
  # Go through these models
  for (g in gs) {

    m <- findModel(model, g)

    z <- z.effects(m, X, Y)
    zs[rownames(z), g] <- z

    # FIXME: sign switch not taken into account in m@W!
    #wx <- m@W$X[rownames(m@W$X) %in% region.genes,]
    #wy <- m@W$Y[rownames(m@W$Y) %in% region.genes,]
    we <- W.effects(m, X, Y)
    wx <- we$X
    wy <- we$Y
    comx <- intersect(names(wx), rownames(wxs))
    comy <- intersect(names(wy), rownames(wys))
    wxs[comx, g] <- wx[comx]
    wys[comy, g] <- wy[comy]
    
  }

  # FIXME: add some check that wxs (and wys) mostly give very similar results
  # otherwise misinterpretations will occur
  #sort(cor(wxs, use = "pairwise.complete.obs")[1,])

  W <- list(X = rowMeans(wxs, na.rm = TRUE), Y = rowMeans(wys, na.rm = TRUE))
  list(z = rowMeans((as.data.frame(zs)), na.rm = TRUE), W = W )
  
}



order.feature.info <- function (feature.info) {

  # Remove genes with no location information from the annotations
  nainds <- is.na(feature.info$chr) | (is.na(feature.info$chr) | is.na(feature.info$arm))
  if (sum(nainds) > 0) {
    feature.info <- feature.info[!nainds, ]
  }
  
  # Order by chromosomal locations
  if ("X" %in% as.character(feature.info$chr)) {
    feature.info$chr[as.character(feature.info$chr) == "X"] <- 23
  }
  if ("Y" %in% as.character(feature.info$chr)) {
    feature.info$chr[as.character(feature.info$chr) == "Y"] <- 24
  }  
  
  feature.info.ordered <- NULL
  chrs <- sort(unique(feature.info$chr))
  arms <- sort(unique(feature.info$arm))
  for (chr in chrs) {
    if (!is.null(arms)) {
      for (arm in arms) {
        #arm.info <- subset(feature.info, chr == chr & arm == arm)
        arm.info <- feature.info[feature.info$chr == chr & feature.info$arm == arm,]
        o <- order(arm.info$loc)
        feature.info.ordered <- rbind(feature.info.ordered, arm.info[o, ])
      }
    } else {
      arm.info <- feature.info[feature.info$chr == chr,]      
      o <- order(arm.info$loc)
      feature.info.ordered <- rbind(feature.info.ordered, arm.info[o, ])
    }
  }
  feature.info.ordered
}



augment.region <- function (region.genes, gene.info) {

  # Check gene.info for genes that are, according to the annotations
  # in gene.info, within the same region than the listed
  # region.genes. Useful for instance when some genes have been left
  # out from modeling due to required one-to-one match between the two
  # data sets.

  # Pick annotations for this region
  #reg.info <- subset(gene.info, geneName %in% region.genes)
  reg.info <- gene.info[gene.info$geneName %in% region.genes, ]

  chrs <- reg.info$chr
  if ( length(unique(chrs)) > 1 ) { stop("Multiple chromosomes listed for the region!") }

  arms <- reg.info$arm
  if ( length(unique(arms)) > 1 ) { stop("Multiple arms listed for the region!") }
  
  reg.start <- as.numeric(min(reg.info$loc))
  reg.end <- as.numeric(max(reg.info$loc))

   # List all genes that are within the same region
  inds <- (gene.info$chr %in% chrs & gene.info$arm %in% arms & gene.info$loc >= reg.start & gene.info$loc <= reg.end )
  rg <-gene.info$geneName[inds]

  # Ensure that also original region.genes are listed
  # even if they are not in gene.info
  union(rg, region.genes)
  
}

