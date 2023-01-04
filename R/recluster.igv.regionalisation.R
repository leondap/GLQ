recluster.igv.regionalisation<-function(fasta,spec,longi,lati, longr=NULL, latr=NULL, minimum=3, minarea=2, nidw=15, gsttab=NULL, GST="all", alpha=0.1, power=2, powerW=1, powerD=1, poweridw=2, method="dist", methodclust="ward.D"){
	res<-NULL
	species<-unique(spec)
	specloc<-aggregate(rep(1,length(longi))~longi+lati+spec,FUN="sum")
	specarea<-aggregate(rep(1,nrow(specloc))~specloc[,3],FUN="sum")

	#Identify species occurring in the minimum number of area and select the data
	speciessel<-specarea[which(specarea[,2]>=minarea),1]
	sel<-which(spec%in%speciessel)
	longi<-longi[sel]
	lati<-lati[sel]
	fasta<-fasta[sel]
	species<-speciessel
	spec<-spec[sel]
	loca<-paste(longi,"_",lati)
	areeun<-unique(loca)
	areetot<-unique(cbind(longi,lati))
	aree<-as.data.frame(cbind(c(1:nrow(areetot)),areetot,areeun))
	loc<-NULL

	for(ar in 1:nrow(aree)){
		wh<-which(loca==areeun[ar])
		loc[wh]<-ar
	}

	if(!(is.null(longr)) & !(is.null(latr))){
		barlong<-NULL
		barlat<-NULL
		longru<-longr[sel]
		latru<-latr[sel]
		barycenters<-cbind(aggregate(longru~loca,FUN="mean"),aggregate(latru~loca,FUN="mean")[,2])
		for (sit in 1:nrow(aree)){
			wh<-which(barycenters[,1]==aree[sit,4])
			barlong[sit]<-barycenters[wh,2]
			barlat[sit]<-barycenters[wh,3]
		}
		aree<-cbind(aree,barlong,barlat)
	}

# Create the dissimilarity matrices for each species (Fig. S1, panel 4)

	matrixGst<-array(NA,dim=c(nrow(aree),nrow(aree),length(species)))
	areatabledissGst<-array(NA,dim=c(nrow(aree),nrow(aree)))
	for (sp in 1:length(species)){
		sel<-which(spec==species[sp])
		fastared<-fasta[sel]
		datamatsp<-loc[sel]
		dismatsp<-as.matrix(dist.dna(fastared, model = "raw",pairwise.deletion = TRUE))
		areamatsp<-unique(datamatsp)

# The contribution of a species can be weighted by its GST on the overall distribution

		if(GST=="all"){
			wei<-1
		}
		if(GST=="gd"){
			if(is.null(gsttab)){
				gstsp<-recluster.fst(as.dist(dismatsp),datamatsp,setzero=T,setnazero=T)$Gst
			}else{
				gstsp<-gsttab[which(names(gsttab)==species[sp])]
			}
			wei<-(1-gstsp)
		}
		if(GST=="st"){
			if(is.null(gsttab)){
				gstsp<-recluster.fst(as.dist(dismatsp),datamatsp,setzero=T,setnazero=T)$Gst
			}else{
				gstsp<-gsttab[which(names(gsttab)==species[sp])]
			}
			wei<-(gstsp)
		}
		if (length(areamatsp)>1){
			for (n in 1: (length(areamatsp)-1)){
				for (m in (n+1):length(areamatsp)){
					dista<-dismatsp[c(which(datamatsp==areamatsp[n]), which(datamatsp==areamatsp[m])),c(which(datamatsp==areamatsp[n]), which(datamatsp==areamatsp[m]))]
					popu<-c(rep(1,length(which(datamatsp==areamatsp[n]))),rep(2,length(which(datamatsp==areamatsp[m]))))
					value<-((recluster.fst(dista,popu,setzero=T,setnazero=T)$Ht)^power)*(wei)
					matrixGst[areamatsp[n], areamatsp[m],sp]<-value
					matrixGst[areamatsp[m], areamatsp[n],sp]<-value
				}
			}
		}
	}

#Compute the average between individual distance matrices (Figure S1, panel 5)

	count<-areatabledissGst
	for (r in 1: nrow(aree)){
		for (c in 1: nrow(aree)){
			count[r,c]<-length(which(!is.na(matrixGst[r,c,])))
			if(count[r,c]>1){
				if(method=="dist"){
					areatabledissGst[r,c]<-(sum(matrixGst[r,c,],na.rm=T)^(1/power))/count[r,c]
					areatabledissGst[c,r]<-(sum(matrixGst[r,c,],na.rm=T)^(1/power))/count[r,c]
				}
				if(method=="mean"){
					areatabledissGst[r,c]<-mean(matrixGst[r,c,],na.rm=T)
					areatabledissGst[c,r]<-mean(matrixGst[r,c,],na.rm=T)
				}
				if(method=="median"){
					areatabledissGst[r,c]<-median(matrixGst[r,c,],na.rm=T)
					areatabledissGst[c,r]<-median(matrixGst[r,c,],na.rm=T)
				}
			}
		}
	}

# Select the areas with enough data for igw (Fig S1, panel 6)

	matnew<-areatabledissGst
	rich<-aggregate(rep(1,length(loc))~loc+spec,FUN="sum")
	richness<-aggregate(rep(1,nrow(rich))~rich[,1],FUN="sum")
	sitiguida<-which(richness[,2]>=minimum)
	guida<-rep(0,nrow(richness))
	guida[sitiguida]<-1
	rich<-richness[,2]
	areeg<-cbind(aree,rich,guida)
	distance<-rep(1,nrow(areeg))
	aree2<-NULL
	k<-1

# Compute igw and fill the empty areas (Fig S1, panel 7)

	for (qua in 1:nrow(matnew)){
		present<-matnew[,qua]
		tabella<-cbind(areeg,present)
		tabella<-tabella[-qua,]
		tabella<-tabella[which(tabella$guida==1),]
		if(length(which(!is.na(tabella$present)))>=nidw){
			tabella<-tabella[complete.cases(tabella),]
			res1 <- idw(as.numeric(tabella$present), cbind(as.numeric(tabella[,2]),as.numeric(tabella[,3])),cbind(as.numeric(areeg[,2]),as.numeric(areeg[,3])),p=poweridw)
			distance<-cbind(distance,res1)
			aree2[k]<-as.character(areeg[qua,1])
			k<-k+1
		}
	}

# Adjust the matrix and make it symmetric (Fig S1, panel 8)
	aree2<-as.numeric(aree2)
	distance<-distance[aree2,]
	distance<-distance[,-1]
	areefin<-sitiguida[which(sitiguida %in% aree2)]
	sel<-which(rownames(distance)%in%areefin)
	dist<-distance[sel,sel]
	mat<-matrix(NA,nrow(dist),nrow(dist))
	for(row in 1:nrow(mat)){
		for(col in 1:ncol(mat)){
			mat[row,col]<-(dist[row,col]+dist[col,row])/2
		}
	}
# Compute geographic distances among areas and compute clustering with soft contiguity constraints  (Fig S1, panel 8)
	dista<-as.dist(mat)
	select<-as.numeric(rownames(dist))
	coord<-areeg[select,]
	if(!(is.null(longr)) & !(is.null(latr))){
		distgeo<-(dist(coord[,5:6]))^powerD
	}else{
		distgeo<-(dist(coord[,2:3]))^powerD
	}
	clust<-hclustgeo2(dista,distgeo,alpha=alpha,w=coord$rich^powerW,method=methodclust)

# Organise the results

	res$richness<-coord$rich
	res$tree <- clust$solution
	res$weightmatrix<-clust$weighteddiss
	res$dist<-dista
	res$distgeo<-distgeo
	res$coord<-cbind(as.numeric(coord[,2]),as.numeric(coord[,3]))
	res$originalcoord<-cbind(longi,lati)
	res$originalsequences<-fasta
	res$originalspec<-spec
	res$minimum<-minimum
	res$minarea<-minarea
	res$gsttab<-gsttab
	res$GST<-GST
	res$alpha<-alpha
	res$power<-power
	res$powerW<-powerW
	res$method<-method
	return(res)
}



hclustgeo2<-function (D0, D1 = NULL, alpha = 0, scale = TRUE, wt = NULL, method="ward.D") 
{
res<-NULL    
if (class(D0) != "dist") 
        stop("DO must be of class dist (use as.dist)", 
            call. = FALSE)
    if (!is.null(D1) && (class(D1) != "dist")) 
        stop("D1 must be of class dist (use as.dist)", 
            call. = FALSE)
    n <- as.integer(attr(D0, "Size"))
    if (is.null(n)) 
        stop("invalid dissimilarities", call. = FALSE)
    if (is.na(n) || n > 65536L) 
        stop("size cannot be NA nor exceed 65536", call. = FALSE)
    if (n < 2) 
        stop("must have n >= 2 objects to cluster", call. = FALSE)
    if (!is.null(D1) && length(D0) != length(D1)) 
        stop("the two dissimilarity structures must have the same size", 
            call. = FALSE)
    if ((max(alpha) > 1) || (max(alpha) < 0)) 
        stop("Values alpha must be in [0,1]", call. = FALSE)
    if ((scale == TRUE) && (!is.null(D1))) {
        D0 <- D0/max(D0)
        D1 <- D1/max(D1)
    }
    delta0 <- wardinit(D0, wt)
    if (!is.null(D1)) 
        delta1 <- wardinit(D1, wt)
    else delta1 <- 0
    delta <- (1 - alpha) * delta0 + alpha * delta1
    res$solution <- hclust(delta, method = method, members = wt)
	res$weighteddiss<-delta
    return(res)
}




