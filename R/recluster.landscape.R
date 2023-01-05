recluster.landscape<-function(mat,units=NULL,dist,transcorr,map_alt,map1,minsea=3){

#Obtain genetic distances (Fig. S4, step 2)

	gendistances<-as.matrix(dist)
	if(!is.null(units)){
		tabuni<-gendistances
		for(r in 1:nrow(gendistances)){
			for(c in 1:nrow(gendistances)){
				if(units[c]==units[r]){
					tabuni[c,r]<-0
				}else{
					tabuni[c,r]<-1
				}
			}
		}
	gendistances<-gendistances*tabuni
	}
#Add a small random position to avoid complete overlap among specimens

	mat[,1]<-mat[,1]+runif(nrow(mat), min = 0.001, max = 0.002)
	mat[,2]<-mat[,2]+runif(nrow(mat), min = 0.001, max = 0.002)
		
#Create the Delaunay triangulation (Fig. S4, step 3)

	try <- deldir(mat[,1],mat[,2])
	maxdat<-max(max(try$delsgs[,5]),max(try$delsgs[,6]))
	if(maxdat<nrow(mat))
		{stop("Not all objects in triangulation consider collinearity")}

#attribute to each edge its values (midpoint, genetic distances, length of the point, maximum and minimum altitude)

	midpointx<-NULL
	midpointy<-NULL
	gendist<-NULL
	lengthpath<-NULL
	lengthsea<-NULL
	seapointx<-NULL
	seapointy<-NULL
	altimin<-NULL
	altimax<-NULL

	for (i in 1:nrow(try$delsgs)){
		gendist[i]<-gendistances[try$delsgs[i,5],try$delsgs[i,6]]
		if(abs(try$delsgs[i,1]-try$delsgs[i,3])<0.043 & abs(try$delsgs[i,2]-try$delsgs[i,4])<0.043){
			try$delsgs[i,1]<-try$delsgs[i,1]+((0.043-abs(try$delsgs[i,1]-try$delsgs[i,3]))*sign(try$delsgs[i,1]-try$delsgs[i,3]))
			try$delsgs[i,2]<-try$delsgs[i,2]+((0.043-abs(try$delsgs[i,2]-try$delsgs[i,4]))*sign(try$delsgs[i,2]-try$delsgs[i,4]))
		}
		patch<-shortestPath(transcorr, c(try$delsgs[i,1],try$delsgs[i,2]),c(try$delsgs[i,3],try$delsgs[i,4]),output="SpatialLines")
		#lines(patch,col="red")
		pa<-as(patch, "sf")
		nume<-round(st_length(pa)*30.36,0)
		punti<-st_sample(pa,nume,type = "regular")
		coords<-(st_coordinates(punti)[,1:2])
		coords<-rbind(c(10,10),coords)
		value<-c(as.vector(extract(map1,coords)))
		altitudes<-as.data.frame(as.vector(extract(map_alt,coords)))
		val<-cbind(coords,value,altitudes)
		altimin[i]<-min(val[,4],na.rm=T)
		altimax[i]<-max(val[,4],na.rm=T)
		lengthpath[i]<-nume
		if(nume==1){
			midpointx[i]<-(try$delsgs[i,1]+try$delsgs[i,3])/2
			midpointy[i]<-(try$delsgs[i,2]+try$delsgs[i,4])/2
			seapointx[i]<-(try$delsgs[i,1]+try$delsgs[i,3])/2
			seapointy[i]<-(try$delsgs[i,2]+try$delsgs[i,4])/2
			lengthsea[i]<-0
		}
		if(nume>1){
			val<-val[-1,]
			tabe<-with(rle(val[,3]), {
 			 	ok <- values==1
 			 	ends <- cumsum(lengths)[ok]
  				starts <- ends - lengths[ok] + 1
  				resi<-as.data.frame(cbind(starts, ends))
				return(resi)
			})
			tabe<-cbind(tabe,(tabe[,2]-tabe[,1]+1))
			tabe<-tabe[which(tabe[,3]>=minsea),]
			line<-which.max(tabe[,3])[1]
			midpointx[i]<-as.numeric(val[round((lengthpath[i]+0.0001)/2,0),1])
			midpointy[i]<-as.numeric(val[round((lengthpath[i]+0.0001)/2,0),2])
			lengthsea[i]<-sum(tabe[,3])
			if(lengthsea[i]>0){
				seapointx[i]<-as.numeric(val[round((tabe[line,1]+tabe[line,2])/2,0),1])
				seapointy[i]<-as.numeric(val[round((tabe[line,1]+tabe[line,2])/2,0),2])
			}
			if(lengthsea[i]==0){
				seapointx[i]<-midpointx[i]
				seapointy[i]<-midpointy[i]
			}
		}
	}
	table<-cbind(try$delsgs,gendist,midpointx,midpointy,lengthpath,seapointx,seapointy,lengthsea,altimin,altimax)
	return(table)
}
