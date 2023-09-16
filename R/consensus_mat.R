
consensus_mat<-function(mat){
	resu<-NULL
	res<-matrix(NA,nrow(mat),nrow(mat))
	for(c in  1:(nrow(mat)-1)){
		for(r in (c+1):(nrow(mat))){
			res[r,c]<- 1-(length(which((mat[c,]-mat[r,])==0))/length(which(!(is.na(mat[c,])) & !(is.na(mat[r,])))))
		}
	}
	resu$dist<-as.dist(res)
	return(resu)
}
