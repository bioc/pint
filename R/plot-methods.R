plot.GeneDependencyModel <- function(x, X, Y = NULL, ann.types = NULL, ann.cols = NULL, legend.x = 0, legend.y = 1, legend.xjust = 0, legend.yjust = 1, order=FALSE, cex.z = 0.6, cex.WX = 0.6, cex.WY = 0.6,...){
                                 
#X <- geneExp; Y <- geneCopyNum; ann.types = NULL; ann.cols = NULL; legend.x = 0; legend.y = 1; legend.xjust = 0; legend.yjust = 1; order=FALSE; cex.z = 0.6; cex.WX = 0.6; cex.WY = 0.6

  model <- x
  
  # Check that both data sets are given for models from 2 data sets
  if (!is.null(getW(model)$X)){
    if (missing(X) || missing(Y)) {
      stop("Original data needed as 'X' and 'Y' arguments")
    }
  } else {   # Check that data set is given for models from 1 data set
    if (missing(X)) {
      stop("Original data needed as 'X' argument")
    }
  }
  z <- z.effects(model, X, Y)[, 1]
  W <- W.effects(model, X, Y)

  if ( order ){
    z.order <- order(z)
    z <- z[as.numeric(z.order)]
  }
  
  # Colors for different annotation types
  cols <- 'grey'
  if (!is.null(ann.types)) {
    if (length(ann.types) != ncol(X$data)) {
      warning("Length of ann.types doesn't match samples in data")
    }
    else {
      labels <- levels(ann.types) 
      types <- max(as.integer(ann.types),na.rm=TRUE)
      if (any(is.na(ann.types))) {
        types <- types + 1
	labels <- c(labels,"NA")
      }
      if (is.null(ann.cols)) {
        ann.cols <- gray(0:types / types)[1:types]
      }
      ann.int <- as.integer(ann.types)
      ann.int <- ifelse(is.na(ann.int),types,ann.int)
      cols <- ann.cols[ann.int]
      if (order) {
        cols <- cols[as.numeric(z.order)]
	ann.int <- ann.int[as.numeric(z.order)]
      }
    }
  }

  def.par <- par(no.readonly = TRUE) # save default

  # Outer margins and layout
  par(oma = c(0.5,0.5,2.5,0.5))

  # For models from 2 data sets
  if (!is.null(Y)){
    if (length(W$X) == 1){
      layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), c(1,3),c(8,6))
    }
    else {
      layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), c(2,2),c(8,6))
    }
  } else { # For models from 1 data set
    layout(matrix(c(1,2),2,1,byrow = TRUE), 1,c(8,6))
  }

  # Inner margins
  par(mar = c(3.1,4.1,3.1,0),cex.main = 1,cex.axis = 1)

  # Barplots for z
  barplot(z,main="Sample effects",col = cols,xlab = "Sample", ylab = "Weight",las = 2,cex.names=cex.z)  

  # Put legend 
  if(!is.null(ann.types)){
    legend(legend.x, legend.y, labels, cex=0.8, ann.cols, xjust=legend.xjust, yjust=legend.yjust);
  }

  par(mar = c(5.1,4.1,3.1,0.5),cex.main = 1,cex.axis = 1)
  # Barplots for Ws
  # For models from 2 data sets
  if (!is.null(Y)){
    if(length(W$X) == 1){
      barplot(t(W$X),main="Variable effects,\n(first data set)",ylab="Weight",las=2,cex.names=cex.WX)
    }
    else {
      barplot(t(W$X),main="Variable effects (first data set)",ylab="Weight",las=2,cex.names=cex.WX)
    }
    barplot(t(W$Y),main="Variable effects (second data set)",ylab="Weight",las=2,cex.names=cex.WY)
  }
  # For models from 1 data set
  else {
    barplot(t(W$total),main="Variable effects",ylab="Weight",las=2,cex.names=cex.WX)
  }

  #Title
  title <- paste(getModelMethod(model),"model around gene",getGeneName(model))
  if(length(getLoc(model)) > 0)
    title <- paste(title,"at",(getLoc(model)/1e6)," ") # was Mbp but removed since user may give locations with kbp or other measure

  mtext(title, NORTH<-3, line=0, adj=0.5, cex=1.2, outer=TRUE)

  par(def.par)
}


plot.ChromosomeModels <-
function(x, hilightGenes = NULL, showDensity = FALSE, showTop = 0, topName = FALSE,
	type = 'l', xlab = 'gene location', ylab = 'dependency score',
	main = NULL,
	pch = 20, cex = 0.75, tpch = 3, tcex = 1, xlim = NA, ylim = NA,...){

  models <- x
  #No printing without any models
  if(isEmpty(models)) 
    return()

  chr <- getChromosome(models)
  arms <- getArm(models)
  scores <- getScore(models)
  locs <- getLoc(models)
  geneNames <- getGeneName(models) 
  if (!any(arms == "")){  
    #scores, locations and gene names for separate arms
    qArm <- models[['q']]
    pArm <- models[['p']]
    qscores <- getScore(qArm)
    qlocs <- getLoc(qArm)
    pscores <- getScore(pArm)
    plocs <- getLoc(pArm)
    #limits for plotting area
    ylim <- c(0,max(c(pscores,qscores)))	
    xlim <- c(min(c(plocs/1e6,qlocs/1e6)),max(c(qlocs/1e6,plocs/1e6)))
  } else {
    #limits for plotting area
    ylim <- c(0,max(scores))	
    xlim <- c(min(locs/1e6),max(locs/1e6))
  }

  # Text for plot
  if (is.null(main)) {
    main <- paste('Dependency score for chromosome ', chr, sep = '')
    if (all(arms == 'p')) main <- paste(main, 'p', sep='')
    if (all(arms == 'q')) main <- paste(main, 'q', sep='')
  }    

  # Plotting dep scores

  if (any(arms == "")){  
    # The whole chromosome
    pl <- plot((locs/1e6), scores, type = type, xlab = xlab, xlim = xlim,
               ylab = ylab, main = main, ylim = ylim, ...)
  } else {
    # Plotting both arms separately

    #p Arm plot
    pl <- plot((plocs/1e6), pscores, type = type, xlab = xlab, xlim = xlim,
               ylab = ylab, main = main, ylim = ylim, ...)
	
    pl <- par(new = TRUE)
	
    #q Arm plot
    pl <- plot((qlocs/1e6), qscores, type = type, xlab = xlab, xlim = xlim,
               ylab = ylab, main = main, ylim = ylim, ...)
	
    pl <- par(new = FALSE)

	}
  if(showTop > 0){
		
    d <- data.frame(scores,locs)
    d <- d[order(scores,decreasing=TRUE),]
		
    #Put vertical dashed line in the middle of showTop highest and showTop+1 highest scores
    limit = (d$score[showTop]+d$score[showTop+1])/2
    abline(h=limit,lty=2)
    topMods <- topModels(models,showTop) 
    #Draw points and add ranking number
    for(i in 1:showTop){
      points(d$locs[i]/1e6,d$scores[i],pch=tpch,cex = tcex) 
      if(topName)
        text(d$locs[i]/1e6,d$scores[i],getGeneName(topMods[[i]]),pos=4,cex=tcex)
      else
        text(d$locs[i]/1e6,d$scores[i],as.character(i),pos=2)
    }
  }
	
  if (!is.null(hilightGenes)){

    matches <- match(hilightGenes, geneNames)
    points(locs[matches]/1e6,scores[matches],pch=pch,cex=cex)
  }
  #lines to bottom to show gene density
  if(showDensity){
    heightCoefficient <- 100
    for(i in 1:length(locs)){
      lines( c((locs[i]/1e6) , (locs[i]/1e6)) , c(0,ylim[2]/heightCoefficient))
    }
  }
}

plot.GenomeModels <-
function(x, hilightGenes = NULL, showDensity = FALSE, showTop = 0, topName = FALSE, onePlot = FALSE, type = 'l',
	 ylab = "Dependency Scores", xlab = "Gene location (chromosome)", 
         main = "Dependency Scores in All Chromosomes",
	 pch = 20, cex = 0.75, tpch = 20, tcex = 0.7, 
	 mfrow = c(5,5), mar = c(3,2.5,1.3,0.5), ps = 5, mgp = c(1.5,0.5,0),ylim=NA,...){

  models <- x
  
  # Save default pars
  op <- par(no.readonly = TRUE)
  
  if (any(is.na(ylim))) {
    ylim = c(0,getScore(topModels(x,1)[[1]]))
  }

  if (!onePlot) {
      if (showTop > 0){
   

      scores <- vector()
      locs <- vector()
      chrs <- vector()
      genes <- vector()
      for (i in 1:24){
        if (isEmpty(x[[i]]))
          next   
        scores <- c(scores, getScore(getPArm(x[[i]])), getScore(getQArm(x[[i]])))
        locs <- c(locs, getLoc(getPArm(x[[i]])), getLoc(getQArm(x[[i]])))
        chrs <- c(chrs, rep(i,getModelNumbers(getPArm(x[[i]])) + getModelNumbers(getQArm(x[[i]]) )))
        genes <- c(genes, getGeneName(getPArm(x[[i]])), getGeneName(getQArm(x[[i]])))
      }

      d <- data.frame(scores,locs,chrs,genes)
      d <- d[order(scores,decreasing=TRUE),]
      # Put vertical dashed line in the middle of showTop highest and showTop+1 highest scores
      limit = (d$scores[showTop]+d$scores[showTop+1])/2
    }
    
    pl <- par(mfrow = mfrow, mar = mar, ps = ps, mgp = mgp)
    for(i in 1:24){
      pl <- plot.ChromosomeModels(models[[i]], hilightGenes, showDensity, ylim=ylim, type=type,
                                  showTop=FALSE,hilightGenes=hilightGenes, showDensity=showDensity,
                                  ylab = ylab, xlab = xlab, main = main, pch = pch, cex = cex, ...)

      if (showTop > 0){
        if(length(which(d$chrs == i)) == 0)
          next
        abline(h=limit,lty=2)
        for (j in which(d$chrs == i)){
          if(j > showTop)
            next
          points(d$locs[j]/1e6,d$scores[j],pch=tpch,cex = tcex)
          if (topName)
            text(d$locs[j]/1e6,d$scores[j],d$genes[j],pos=4,cex=tcex)
          else
            text(d$locs[j]/1e6,d$scores[j],as.character(j),pos=2,cex=tcex)
        }
      }
    }
  }

  else {     
    plocs <- list()
    pscores <- list()
    qlocs <- list()
    qscores <- list()
    drawLoc <- vector()
    locStart <- 0
    drawScore <- vector()
    geneNames <- vector()
    tickMarks <- vector()
    tickWidth <- 50*1e6
    tickLocs <- vector()
    lineLocs <- vector()

    for (i in 1:24){
      if (isEmpty(x[[i]]))
        next
      plocs[[i]] <- getLoc(getPArm(x[[i]]))
      pscores[[i]] <- getScore(getPArm(x[[i]]))
      qlocs[[i]] <- getLoc(getQArm(x[[i]]))
      qscores[[i]] <- getScore(getQArm(x[[i]]))
      drawLoc <- c(drawLoc, locStart+plocs[[i]], locStart+qlocs[[i]])
      drawScore <- c(drawScore, pscores[[i]], qscores[[i]])
      geneNames <- c(geneNames, getGeneName(getPArm(x[[i]])), getGeneName(getQArm(x[[i]])))
      lineLocs <- c(lineLocs, (locStart)/1e6)
      tickLocs <- c(tickLocs, (locStart+max(drawLoc))/2e6)
      tickMarks <- c(tickMarks, i)
      locStart <- locStart + max(plocs[[i]], qlocs[[i]])
    }
    lineLocs <- c(lineLocs, (locStart)/1e6)

    plot(0,0,type = type, xlab=xlab, xlim=c(0,max(drawLoc/1e6)),
        ylab=ylab, main=main, ylim = ylim, xaxt = 'n',...)

    locStart <- 0
    for(i in 1:24){
      if (isEmpty(x[[i]]))
        next
      lines((locStart+ plocs[[i]])/1e6, pscores[[i]], type=type) 
      lines((locStart+ qlocs[[i]])/1e6, qscores[[i]], type=type) 
      locStart <- locStart + max(plocs[[i]], qlocs[[i]])
      abline(v=lineLocs[i+1])
    }
    #plot((drawLoc/1e6), drawScore, type = type, xlab=xlab,
    #    ylab=ylab, main=main, ylim = ylim, xaxt = 'n',...)

    axis(1,at=tickLocs,labels=tickMarks,cex.axis=0.7,tick=FALSE)
  
    if(showTop > 0){
   
      d <- data.frame(drawScore,drawLoc)
      d <- d[order(drawScore,decreasing=TRUE),]

      # Put vertical dashed line in the middle of showTop highest and showTop+1 highest scores
      limit = (d$drawScore[showTop]+d$drawScore[showTop+1])/2
      abline(h=limit,lty=2)
      topMods <- topModels(models,showTop) 
      # Draw points and add ranking number
      for(i in 1:showTop){
        points(d$drawLoc[i]/1e6,d$drawScore[i],pch=tpch,cex = tcex)
	  if(topName)
	    text(d$drawLoc[i]/1e6,d$drawScore[i],getGeneName(topMods[[i]]),pos=4,cex=tcex)
	  else
          text(d$drawLoc[i]/1e6,d$drawScore[i],as.character(i),pos=2,cex=tcex)
      }
    }

    if (!is.null(hilightGenes)){
      # indices of cancer genes
      matches <- match(hilightGenes, geneNames)
      # print(matches)
      points(drawLoc[matches]/1e6,drawScore[matches],pch=pch,cex=cex)
    }

    # lines to bottom to show gene density
    if (showDensity){
      heightCoefficient <- 100
      for(i in 1:length(drawLoc)){
        lines( c((drawLoc[i]/1e6) , (drawLoc[i]/1e6)) , c(0,ylim[2]/heightCoefficient))
      }
    }
  }
  # Reset default pars
  par(op)
}


     




