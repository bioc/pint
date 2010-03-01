matrix.sqrt <-
function (A) {
	#Solve square root (A^(1/2)) of a diagonalizable nxn-matrix A
	#First diagonalize A, then compute sqrt for the diagonal elements
	#of the diagonalized matrix Adot. Then compute matrix sqrt from these
	#Eigenvectors of A

	#Author: Leo Lahti
	#Date: 13.2.2007
	#Comment: I did not easily find suitable function in R although there might well be a better optimized one. 

	e<-eigen(A) 
	V<-e$vector
	#try inverse the eigvec matrix. 
	#if not diagonalizable then some error will occur from try I think
	Vinv<-try(solve(V)) 
	D <- diag(e$values)

	V%*%sqrt(D)%*%Vinv
}

