plot.DependencyModel <- function(x, X, Y, ...){

  model <- x
  if(missing(X) || missing(Y)) {
    stop("Original data needed")
  }
  z = z.projection(model,X,Y)
  W = W.projection(model,X,Y)

  # Outer margins for title of dependency model
  par(oma = c(0.5,0.5,2.5,0.5))
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  # Inner margins
  par(mar = c(5.1,4.1,3.1,0.5),cex.main = 1,cex.axis = 1)
  barplot(t(z),main="Sample effects",xlab = "Sample", ylab = "Weight",las = 2,cex.names=0.5)
  barplot(t(W$X),main="Variable effects (first data set)",ylab="Weight",las=2,cex.names=0.5)
  barplot(t(W$Y),main="Variable effects (second data set)",ylab="Weight",las=2,cex.names=0.5)

  #Title
  title = paste(getModelMethod(model),"model around gene",getGeneName(model))
  if(length(getLoc(model)) > 0)
    title = paste(title,"at",(getLoc(model)/1e6),"Mbp")

  mtext(title, NORTH<-3, line=0, adj=0.5, cex=1.2, outer=TRUE)
}

plot.ChromosomeArmModels <-
function (x, hilightGenes = NULL, showDensity = FALSE, showTop = 0,
	type = 'l', xlab = 'gene location (Mbp)', ylab = 'dependency score',
	main = paste('Dependency score for chromosome ', chr, arm, sep = ''),
	pch = 20, cex = 0.75, tpch = 3, tcex = 1, ylim = NA, ...){
	models <- x
	
	#No printing without any models
	if(getModelNumbers(models) == 0) return()
	
	scores <- getScore(models)
	locs <- getLoc(models)
	chr <- getChromosome(models)
	arm <- getArm(models)
	geneNames <- getGeneName(models)

	if(all(is.na(ylim))){
		ylim <- c(0,max(scores))
	}

	plot((locs/1e6), scores, type = type, xlab=xlab, 
		ylab=ylab, main=main, ylim = ylim,...)
	
	if(showTop > 0){
		topModels <- findHighestModels(models,showTop)
		
		d <- data.frame(scores,locs)
		d <- d[order(scores,decreasing=TRUE),]
		
		#Put vertical dashed line in the middle of showTop highest and showTop+1 highest scores
		limit = (d$score[showTop]+d$score[showTop+1])/2
		abline(h=limit,lty=2)
		#Draw points and add ranking number
		for(i in 1:showTop){
			points(d$locs[i]/1e6,d$scores[i],pch=tpch,cex = tcex) 
			text(d$locs[i]/1e6,d$scores[i],as.character(i),pos=2)
		}
	}
	
	if (!is.null(hilightGenes)){

		#indices of cancer genes
		matches <- match(hilightGenes, geneNames)
		#print(matches)
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


plot.ChromosomeModels <-
function(x, hilightGenes = NULL, showDensity = FALSE, showTop = 0,
	type = 'l', xlab = 'gene location (Mbp)', ylab = 'dependency score',
	main = paste('Dependency score for chromosome ', chr, sep = ''),
	pch = 20, cex = 0.75, tpch = 3, tcex = 1, xlim = NA, ylim = NA,...){

	models <- x
	#No printing without any models
	if(getModelNumbers(getPArm(models)) == 0 && getModelNumbers(getQArm(models)) == 0) 
		return()


	pArm <- getPArm(models)
	qArm <- getQArm(models)
	chr <- getChromosome(pArm)
	
	pscores <- getScore(pArm)
	plocs <- getLoc(pArm)
	#p Arm to negative side
	#plocs = -plocs
	pgeneNames <- getGeneName(pArm)
	
	qscores <- getScore(qArm)
	qlocs <- getLoc(qArm)
	qgeneNames <- getGeneName(qArm)

	if(length(qscores) == 0){
		qscores <- 0
		qlocs <- 0
	}
	if(length(pscores) == 0){
		pscores <- 0
		plocs <- 0
	}

	if(all(is.na(ylim))){
		ylim <- c(0,max(max(pscores),max(qscores),1.0))	
	}
	if(all(is.na(xlim))){
		xlim <- c(min(plocs/1e6),max(qlocs/1e6))
	}
	
	#p Arm plot
	pl <- plot((plocs/1e6), pscores, type = type, xlab = xlab, xlim = xlim,
		ylab = ylab, main = main, ylim = ylim, ...)
	
	pl <- par(new = TRUE)
	
	#q Arm plot
	pl <- plot((qlocs/1e6), qscores, type = type, xlab = xlab, xlim = xlim,
		ylab = ylab, main = main, ylim = ylim, ...)
	
	pl <- par(new = FALSE)
	
	if(showTop > 0){
		topModels <- findHighestModels(models,showTop)
		
		scores=c(pscores,qscores)
		locs=c(plocs,qlocs)
		d <- data.frame(scores,locs)
		d <- d[order(scores,decreasing=TRUE),]
		
		#Put vertical dashed line in the middle of showTop highest and showTop+1 highest scores
		limit = (d$score[showTop]+d$score[showTop+1])/2
		abline(h=limit,lty=2)
		#Draw points and add ranking number
		for(i in 1:showTop){
			points(d$locs[i]/1e6,d$scores[i],pch=tpch,cex = tcex) 
			text(d$locs[i]/1e6,d$scores[i],as.character(i),pos=2)
		}
	}
	
	if (!is.null(hilightGenes)){

		#indices of cancer genes in p arm
		pMatches <- match(hilightGenes, pgeneNames)
		#print(pMatches)
		points(plocs[pMatches]/1e6,pscores[pMatches],pch=pch,cex=cex)

		#indices of cancer genes in q arm
		qMatches <- match(hilightGenes, qgeneNames)
		#print(qMatches)
		points(qlocs[qMatches]/1e6,qscores[qMatches],pch=pch,cex=cex)

	}
	#lines to bottom to show gene density
	if(showDensity){
		heightCoefficient <- 100
		for(i in 1:length(plocs)){
			lines( c((plocs[i]/1e6) , (plocs[i]/1e6)) , c(0,ylim[2]/heightCoefficient))
		}
		for(i in 1:length(qlocs)){
			lines( c((qlocs[i]/1e6) , (qlocs[i]/1e6)) , c(0,ylim[2]/heightCoefficient))
		}
	}

	
}

plot.GenomeModels <-
function(x, hilightGenes = NULL, showDensity = FALSE, 
	mfrow = c(5,5), mar = c(3,2.5,1.3,0.5), ps = 5, mgp = c(1.5,0.5,0),...){

	models <- x
	pl <- par(mfrow = mfrow, mar = mar, ps = ps, mgp = mgp)
	
	for(i in 1:24){
		pl <- plot.ChromosomeModels(models[[i]], hilightGenes, showDensity, ...)
	}

}



