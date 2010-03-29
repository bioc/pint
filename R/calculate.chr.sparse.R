calculate.chr.sparse <- function(X, Y, windowSize, chromosome, method = "pSimCCA", params = list()){
	pArm <- calculate.arm.sparse(X, Y, windowSize, chromosome, 'p', method, params)
	qArm <- calculate.arm.sparse(X, Y, windowSize, chromosome, 'q', method, params)
	chromosome <- factor(chromosome, levels = levels(X$info$chr))
	return(new("ChromosomeModels", pArmModels = pArm, qArmModels = qArm, chromosome = chromosome, method = method, params = params))
}