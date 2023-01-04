
select.maximum<-function(datasel,dist,size=0.045,max=3){
res<-NULL
lat<-(floor(datasel[,3]/size)*size)+(size/2)
long<-(floor(datasel[,4]/size)*size)+(size/2)
val<-paste(lat,long)

quant1<-aggregate(rep(1,nrow(datasel))~ val, FUN="sum")
	datasel2<-datasel[1,]
	for(loc in 1 : nrow(quant1)){
		#loc<-1
		which<-which(val==quant1[loc,1])
		if(quant1[loc,2]<(max+1)){
			datasel2<-rbind(datasel2,datasel[which,])
		}
		if(quant1[loc,2]>max){
			dataselt<-datasel[which,]
			sam<-sample(c(1:length(which)))
			dataselt<-dataselt[sam,]
			datasel2<-rbind(datasel2,dataselt[1:max,])
					}
	}


res$data<-datasel2[-1,]	
matchseqs<-match(res$data[,1],datasel[,1])
dista<-as.matrix(dist)[matchseqs,matchseqs]
res$dist<-as.dist(dista)
return(res)
}
