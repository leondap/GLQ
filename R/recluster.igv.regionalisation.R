
recluster.igv.regionalisation<-function(fasta,spec,longi,lati, longr=NULL, latr=NULL, minimum=2, minarea=2, fsttype="Hinter", modeldist="raw", spec_norm=1, gsttab=NULL, GST="all", power=2, method="dist", minshared=3, methodclust="ward.D"){
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
	if(!(is.null(gsttab))){
      m<-match(species,names(gsttab))
	gsttab<-gsttab[m]
	}
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
	matrixnum<-array(NA,dim=c(nrow(aree),nrow(aree),length(species)))
	matrixgsttab<-array(NA,dim=c(nrow(aree),nrow(aree),length(species)))
	matrixfreq<-array(NA,dim=c(nrow(aree),nrow(aree),length(species)))


	areatabledissGst<-array(NA,dim=c(nrow(aree),nrow(aree)))
	for (sp in 1:length(species)){
		sel<-which(spec==species[sp])
		fastared<-fasta[sel]
		datamatsp<-loc[sel]
		dismatsp<-as.matrix(dist.dna(fastared, model = modeldist,pairwise.deletion = TRUE))
		areamatsp<-unique(datamatsp)

# The contribution of a species can be weighted by its GST on the overall distribution

		if(GST=="all"){
			wei<-1
		}
		if(GST=="gd"){
			if(is.null(gsttab)){
				gsttab[sp]<-recluster.fst(as.dist(dismatsp),datamatsp,setzero=T,setnazero=T)$Gst
				gstsp<-gsttab[sp]
			}else{
				gstsp<-gsttab[which(names(gsttab)==species[sp])]
			}
			wei<-(1-gstsp)
		}
		if(GST=="st" | method=="gstst"){
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
					fst<-recluster.fst(dist=dista,vect=popu,setzero=T,setnazero=T)
					if(fsttype=="Hinter") {value<-(fst$Hinter)}
					if(fsttype=="Dst") {value<-(fst$Dst)}
					if(fsttype=="Gst") {value<-(fst$Gst)}
					if(fsttype=="Ht") {value<-(fst$Ht)}
					specinum<-length(popu)
					if(specinum > 15){
						specinum<-16
					}
					matrixGst[areamatsp[n], areamatsp[m],sp]<-value
					matrixGst[areamatsp[m], areamatsp[n],sp]<-value
					matrixnum[areamatsp[n], areamatsp[m],sp]<-specinum
					matrixnum[areamatsp[m], areamatsp[n],sp]<-specinum
					matrixfreq[areamatsp[m], areamatsp[n],sp]<-length(areamatsp)
					matrixfreq[areamatsp[n], areamatsp[m],sp]<-length(areamatsp)

					if(!(is.null(gsttab))){
						matrixgsttab[areamatsp[n], areamatsp[m],sp]<-gsttab[sp]
						matrixgsttab[areamatsp[m], areamatsp[n],sp]<-gsttab[sp]
					}
				}
			}
		}
	}

#Compute the average between individual distance matrices (Figure S1, panel 5)

	count<-areatabledissGst
	GSTvalues<-areatabledissGst

	for (r in 1: nrow(aree)){
		for (c in 1: nrow(aree)){
#r<-1
#c<-46
			specmat<-which(!is.na(matrixGst[r,c,]))
			count[r,c]<-length(specmat)
			GSTvalues[r,c]<-sum(gsttab[specmat])
			denominatore<-count[r,c]
			if(method=="mean"){
			denominatore<-count[r,c]
			}
			
			if(count[r,c]>= minshared){
				if(method=="mean"){
					areatabledissGst[r,c]<-(sum(matrixGst[r,c,]^power,na.rm=T)/denominatore)^(1/power)
					areatabledissGst[c,r]<-(sum(matrixGst[r,c,]^power,na.rm=T)/denominatore)^(1/power)
				}
				if(method=="dist2"){
					areatabledissGst[r,c]<-(sum(matrixGst[r,c,]^power*(log(matrixnum[r,c,],2))^power, na.rm=T)/sum(log(matrixnum[r,c,],2)^power, na.rm=T))^(1/power)
					areatabledissGst[c,r]<-(sum(matrixGst[r,c,]^power*(log(matrixnum[r,c,],2))^power, na.rm=T)/sum(log(matrixnum[r,c,],2)^power, na.rm=T))^(1/power)
				}
				if(method=="gstst"){
					areatabledissGst[r,c]<-(sum((matrixGst[r,c,]^power)*(matrixgsttab[r,c,]^power), na.rm=T)/sum((matrixgsttab[r,c,]^power), na.rm=T))^(1/power)
					areatabledissGst[c,r]<-(sum((matrixGst[r,c,]^power)*(matrixgsttab[r,c,]^power), na.rm=T)/sum((matrixgsttab[r,c,]^power), na.rm=T))^(1/power)

				}
				if(method=="freq"){
					areatabledissGst[r,c]<-(sum((matrixGst[r,c,]^power)*(log(matrixfreq[r,c,])^power), na.rm=T)/sum(log(matrixfreq[r,c,])^power, na.rm=T))^(1/power)
					areatabledissGst[c,r]<-(sum((matrixGst[r,c,]^power)*(log(matrixfreq[r,c,])^power), na.rm=T)/sum(log(matrixfreq[r,c,])^power, na.rm=T))^(1/power)

				}

						
				if(method=="median"){
					areatabledissGst[r,c]<-median(matrixGst[r,c,],na.rm=T)
					areatabledissGst[c,r]<-median(matrixGst[r,c,],na.rm=T)
				}
			}
		}
	}


# Select the areas with complete data

matnew<-areatabledissGst

rich<-aggregate(rep(1,length(loc))~loc+spec,FUN="sum")
richness<-aggregate(rep(1,nrow(rich))~rich[,1],FUN="sum")
sitiguida<-which(richness[,2]>= minimum)
guida<-rep(0,nrow(richness))
guida[sitiguida]<-1
rich<-richness[,2]
areeg<-cbind(aree,rich,guida)
distance<-rep(1,nrow(areeg))

selarea<-which(areeg[,8]==1)
areegs<-areeg[selarea,]
tabgstsel<-areatabledissGst[selarea,selarea]
diag(tabgstsel)<-0
matnewna<-tabgstsel
matnewna[!(is.na(matnewna))]<-0
matnewna[(is.na(matnewna))]<-1
nasum<-rowSums(matnewna)
perc1<-(nasum/length(nasum))*100

selecta<-c(1:nrow(tabgstsel))
righe<-nrow(tabgstsel)
tabgstselu<-tabgstsel

for(add in 1:righe){
matnewna<-tabgstselu
matnewna[!(is.na(matnewna))]<-0
matnewna[(is.na(matnewna))]<-1
nasum<-rowSums(matnewna)
perc<-(nasum/length(nasum))*100
if(sum(perc)==0){
break
}
togli<-which.max(perc)
selecta<-selecta[-togli]
tabgstselu<-tabgstselu[-togli,-togli]
}

areegs2<-areegs[selecta,]
tabgstsel2<-tabgstsel[selecta,selecta]
areegs2<-areegs[selecta,]
removed<-rep(0,nrow(areegs))
removed[selecta]<-1
imputed<-as.dist(tabgstsel2)
	res$dist<-imputed
	res$allareas<-cbind(areegs,removed,perc1)
	res$alldist<-tabgstsel
	res$areas<-areegs2
	res$originalcoord<-cbind(longi,lati)
	res$originalsequences<-fasta
	res$originalspec<-spec
	res$minimum<-minimum
	res$minarea<-minarea
	res$gsttab<-gsttab
	res$GST<-GST
	res$power<-power
	res$method<-method
	return(res)
}


