
recluster.latitudinal.bands<-function(location, identification, sequences, coordinates, populations, range=c(0,90), GST=NULL, partition=c("ALL","GD","ST"), correct=F, width=3,asympt=FALSE, mini_sp=5, iter=10){

#Compute the number of bands	
	interval<-width/2
	start<-range[1]
	slices<-round(((range[2]-range[1])/interval)+1,0)
#Create the empty objects
	res<-NULL
	matrixsp<-array(NA, dim=c(100000,14,iter))
	dataset<-cbind(location,identification,coordinates,populations)
	species<-unique(identification)
	matricerichness<-array(0,dim=c(length(species),slices+1,iter))
	rownames(matricerichness)<-species
	matricehaplotypes<-array(0,dim=c(length(species),slices+1,iter))
	rownames(matricehaplotypes)<-species
	matricehaplotypeslow<-array(0,dim=c(length(species),slices+1,iter))
	rownames(matricehaplotypeslow)<-species
	matricedistances<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matricedistances)<-species
	matriceGD<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matriceGD)<-species
	matriceGDred<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matriceGD)<-species
	matriceST<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matriceGD)<-species
	matriceSTred<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matriceGD)<-species
	matriceALL<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matriceGD)<-species
	matriceALLred<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matriceGD)<-species
	matricedistancesGD<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matricedistancesGD)<-species
	matricedistancesST<-array(NA,dim=c(length(species),slices+1,iter))
	rownames(matricedistancesST)<-species

# Launch the iteration (needed because a subsample of species per band is used to include a similar number of specimens per band)

	for (it in 1:iter){
		for (cut in 0:slices){
			lati<-start+interval*cut
			stripe<-which(dataset[,1]>=lati-interval & dataset[,1]<lati+interval)
			if(length(stripe)>0){
				datalat<-dataset[stripe,]
				speciessp<-aggregate(rep(1,nrow(datalat))~ datalat[,2], FUN="sum")
				for(k in 1:nrow(speciessp)){
					spe<-which(rownames(matricerichness)==speciessp[k,1])
					matricerichness[spe,cut+1,it]<-speciessp[k,2]
				}
			}
		}
	line<-1

#Analyse each species separately

	for (spe in 1:length(species)){

# Identify bands where each species has minimum requirements

		dataspeq<-dataset[which(dataset[,2]==species[spe]),]
		seqspeq<-sequences[which(dataset[,2]==species[spe])]
		stripes<-which(matricerichness[spe,,it]>=mini_sp)
		if(length(stripes)>0){
			minimum<-min(matricerichness[spe,,it][stripes])
			max<-max(matricerichness[spe,,it])
			for (strip in 1:length(stripes)){

#Extract the specimens of this species from each band

				lati<-start+interval*(stripes[strip]-1)
				stripe<-which(dataspeq[,1]>=lati-interval & dataspeq[,1]<lati+interval)
				seqslat<-seqspeq[stripe]
				datalat<-dataspeq[stripe,]

# Calculate the barycentre in each band

				baricenter<-c(mean(datalat[,3]),mean(datalat[,4]))
				coordin<-rbind(baricenter,datalat[,3:4])
				qualiseqs<-sample(1:length(seqslat))[1:minimum]
				baricenter<-c(mean(datalat[,3]),mean(datalat[,4]))
				coordin<-rbind(baricenter,datalat[,3:4])
				molt<-1
				same<-as.vector(dist(datalat$populations))
				same[which(same>0)]<-1
				regi<-unique(datalat$populations)
				regi<-regi[!is.na(regi)]
				GDdis<-NULL
				numero<-0
				for(re in 1:length(regi)){
					piglia<-datalat[which(datalat$populations==regi[re]),]
					if(nrow(piglia)>1){
						numero<-numero+nrow(piglia)
						baricenterGD<-c(mean(piglia[,3]),mean(piglia[,4]))
						coordinGD<-rbind(baricenterGD,piglia[,3:4])
						distanzGD<-as.matrix(earth.dist(coordinGD, dist = TRUE))[1,]
						distanzGD<-distanzGD[2:length(distanzGD)]
						GDdis<-c(GDdis,distanzGD)
					}
				}

				if(!is.null(GST)){
					gst<-GST[which(names(GST)==species[spe])]
					if(partition=="ALL"){
						molt<-1
					}
					if(partition=="GD"){
						molt<-1-gst
					}
					if(partition=="ST"){
						molt<-gst
					}
				}
				tutte<-dist.dna(seqslat, model = "raw",pairwise.deletion = TRUE)
				tuttered<-(dist.dna(seqslat[qualiseqs], model = "raw",pairwise.deletion = TRUE))
				matriceALL[spe,stripes[strip],it]<-(mean(tutte))*molt
				matriceALLred[spe,stripes[strip],it]<-(mean(tuttered))*molt
				if (correct==T){
					matriceALL[spe,stripes[strip],it]<-sqrt((tutte)*(length(seqslat)/(length(seqslat)-1)))*molt
					matriceALLred[spe,stripes[strip],it]<-(mean(tuttered)*(length(qualiseqs)/(length(qualiseqs)-1)))*molt
				}
				intra<-which(same==0)
				matriceGD[spe,stripes[strip],it]<-(mean(as.vector(tutte)[intra]))*molt
				if (correct==T){
					matriceGD[spe,stripes[strip],it]<-(mean(as.vector(tutte)[intra])*(length(qualiseqs)/(length(qualiseqs)-1)))*molt
				}
				distanz<-as.matrix(earth.dist(coordin, dist = TRUE))[1,]
				distanz<-distanz[2:length(distanz)]
				matricedistancesGD[spe,stripes[strip],it]<-sqrt((sum(GDdis^2)/numero))
				inter<-which(same==1)
				matriceST[spe,stripes[strip],it]<-(mean(as.vector(tutte)[inter])*(length(qualiseqs)/(length(qualiseqs)-1)))*molt
				if (correct==T){
					matriceST[spe,stripes[strip],it]<-(mean(as.vector(tutte)[inter])*(length(qualiseqs)/(length(qualiseqs)-1)))*molt
				}
				matricehaplotypeslow[spe,stripes[strip],it]<-length(recluster.haplotypes(seqslat[qualiseqs])$frequency)*molt
				matricedistances[spe,stripes[strip],it]<-sqrt(mean(distanz^2))
				matrixsp[line,1,it]<-species[spe]
				matrixsp[line,2,it]<-matricedistances[spe,stripes[strip],it]
				matrixsp[line,3,it]<-matriceALL[spe,stripes[strip],it]
				matrixsp[line,4,it]<-matricehaplotypeslow[spe,stripes[strip],it]
				matrixsp[line,11,it]<-lati
				matrixsp[line,7,it]<-length(seqslat)
				matrixsp[line,8,it]<-length(qualiseqs)
				matrixsp[line,9,it]<-baricenter[1]
				matrixsp[line,10,it]<-baricenter[2]
				matrixsp[line,6,it]<-matriceALLred[spe,stripes[strip],it]
				matrixsp[line,12,it]<-matricedistancesGD[spe,stripes[strip],it]
				matrixsp[line,13,it]<-matriceGD[spe,stripes[strip],it]
				matrixsp[line,14,it]<-matriceST[spe,stripes[strip],it]
				if(asympt){
					if(length(freq)>1){
					accu<-iNEXT(freq, q=0, datatype="abundance", size=NULL, endpoint=max, knots=40, se=TRUE, conf=0.95, nboot=50)
					matricehaplotypes[spe,stripes[strip],it]<-accu$AsyEst[1,2]*molt
					matrixsp[line,5]<-matricehaplotypes[spe,stripes[strip]]
				}
				if(length(freq)==1){
					matricehaplotypes[spe,stripes[strip],it]<-1*molt
					matrixsp[line,5,it]<-matricehaplotypes[spe,stripes[strip],it]
				}
			}
			line<-line+1
		}
	}
	}
	}

	diversity<-matrix(NA,iter,ncol(matricehaplotypes))
	speciesnum<-diversity
	latitude<-diversity
	distances<-diversity
	reducedhap<-diversity
	quanteslice<-diversity
	reducedGD<-diversity
	reducedhapsd<-diversity
	reducedGDsd<-diversity
	bar_long<-diversity
	bar_lat<-diversity
	centroid<-diversity
	reducedALL<-diversity
	reducedALLred<-diversity
	reducedALLredsd<-diversity
	reducedALLsd<-diversity
	reducedST<-diversity
	reducedSTsd<-diversity
	distancesGD<-diversity
	quanteslice<-NULL
	for (slice in 1:ncol(matricehaplotypeslow)){
		quali<-which(matricehaplotypeslow[,slice,1]>0)
		quanteslice[slice]<-length(quali)
		if(quanteslice[slice]>=mini_sp){
			for(it in 1:iter){
				diversity[it,slice]<-mean(matricehaplotypes[quali,slice,it],na.rm=T)
				speciesnum[it,slice]<-length(quali)
				latitude[it,slice]<-start+interval*(slice-1)
				distances[it,slice]<-mean(matricedistances[quali,slice,it],na.rm=T)
				reducedhap[it,slice]<-mean(matricehaplotypeslow[quali,slice,it],na.rm=T)
				reducedGD[it,slice]<-mean(matriceGD[quali,slice,it],na.rm=T)
				reducedhapsd[it,slice]<-sd(matricehaplotypeslow[quali,slice,it],na.rm=T)
				reducedGDsd[it,slice]<-sd(matriceGD[quali,slice,it],na.rm=T)
				reducedALL[it,slice]<-mean(matriceALL[quali,slice,it],na.rm=T)
				reducedALLsd[it,slice]<-sd(matriceALL[quali,slice,it],na.rm=T)
				reducedALLred[it,slice]<-mean(matriceALLred[quali,slice,it],na.rm=T)
				reducedALLredsd[it,slice]<-sd(matriceALLred[quali,slice,it],na.rm=T)
				reducedST[it,slice]<-mean(matriceST[quali,slice,it],na.rm=T)
				reducedSTsd[it,slice]<-sd(matriceST[quali,slice,it],na.rm=T)
				distancesGD[it,slice]<-mean(matricedistancesGD[quali,slice,it],na.rm=T)
			}
		}
	}
	colnames(matrixsp)<-c("species","Distances","IGV","haplotypes_red","haplotype_asy","IGV_red","Specimens","Specimens_red",
	"Bar_long","Bar_lat","slice","Distances_GD","GD","ST")
	rowgood<-which(complete.cases(matrixsp[,1:2,1]))
	dista<-as.numeric(matrixsp[rowgood,2,1])
	haplotypes_red<-as.numeric(matrixsp[rowgood,4,1])
	IGV_red<-as.numeric(matrixsp[rowgood,6,1])
	Specimens_red<-as.numeric(matrixsp[rowgood,8,1])
	Bar_long<-as.numeric(matrixsp[rowgood,9,1])
	Bar_lat<-as.numeric(matrixsp[rowgood,10,1])
	for(it in 2:iter){
		dista<-cbind(dista,as.numeric(matrixsp[rowgood,2,it]))
		haplotypes_red<-cbind(haplotypes_red,as.numeric(matrixsp[rowgood,4,it]))
		IGV_red<-cbind(IGV_red,as.numeric(matrixsp[rowgood,6,it]))
		Specimens_red<-cbind(Specimens_red,as.numeric(matrixsp[rowgood,8,it]))
		Bar_long<-cbind(Bar_long,as.numeric(matrixsp[rowgood,9,it]))
		Bar_lat<-cbind(Bar_lat,as.numeric(matrixsp[rowgood,10,it]))
	}
	res$separate<-as.data.frame(cbind(matrixsp[rowgood,1,1],rowMeans(dista),matrixsp[rowgood,3,1],rowMeans(haplotypes_red),matrixsp[rowgood,5,1],rowMeans(IGV_red),matrixsp[rowgood,7,1],matrixsp[rowgood,8,1],
	rowMeans(Bar_long),rowMeans(Bar_lat),matrixsp[rowgood,11,1],matrixsp[rowgood,12,1],matrixsp[rowgood,13,1],matrixsp[rowgood,14,1]))
	colnames(res$separate)<-c("species","distances","IGV","haplotypes_red","haplotype_asy","IGV_red","Specimens","Specimens_red",
"Bar_long","Bar_lat","slice","Distances_GD","GD","ST")
	final<-as.data.frame(cbind(colMeans(diversity), colMeans(reducedALL), colMeans(reducedALLsd), colMeans(reducedALLred), colMeans(reducedALLredsd), 
colMeans(reducedhap), colMeans(reducedhapsd),colMeans(latitude),colMeans(speciesnum),colMeans(distances),colMeans(reducedGD),colMeans(reducedGDsd) ,
colMeans(reducedST),colMeans(reducedSTsd),colMeans(distancesGD)))
	res$aggregate<-final[complete.cases(final[,1:2]),]
	return(res)
}



