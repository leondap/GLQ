endemic.mutations<-function(hn,dis,labels,alias=NULL){
	dis<-as.matrix(d)
	link_counts <-c(hn[,1],hn[,2])
	links<-aggregate(rep(1,length(link_counts))~link_counts,FUN="sum")
	connections<-cbind(hn[,1],hn[,2],hn[,3])
	nameshapl<-attr(hn, "labels")
	newvect<-NULL
	if(!(is.null(alias))){
		names<-NULL
		namesd<-NULL
		for(nam in 1:length(nameshapl)){
			nome<-which(alias[,2]==nameshapl[nam])
			nome2<-which(alias[,2]==rownames(dis)[nam])
			names[nam]<-alias[nome,1]
			namesd[nam]<-alias[nome2,1]

		}
	}else{
		names<-nameshapl
		namesd<-rownames(dis)
	}
	rownames(dis)<-namesd
	colnames(dis)<-namesd

	connames<-connections
	for(n in 1:length(names)){
		#n<-1
		w1<-which(connections[,1]==n)
		if(length(w1)>0){
			connames[w1,1]<-names[n]
		}
		w2<-which(connections[,2]==n)
		if(length(w2)>0){
			connames[w2,2]<-names[n]
		}
	}
	connames<-as.data.frame(connames)
	links[,1]<-names
	vect<-rep(NA,length(labels))
	remain<-which(is.na(vect))
	interc<-connames
	link_counts2 <-link_counts
	links2<-links
	names2<-names
	labels2<-labels


	for (giro in 1:1000){
		quali<-which(links2[,2]==1)
		involved<-links2[quali,1]
		selected<-involved[which(involved %in% labels)]
		if(length(selected)>0){
			le<-as.data.frame(cbind(c(interc[,1],interc[,2]),c(interc[,3],interc[,3])))
			if(nrow(le)==2){
				vect[which(labels %in% le[,1])]<-le[1,2]
			}
			if(nrow(le)>2){
				for(g in 1:length(selected)){
					#g<-2
					poni<-which(labels2==selected[g])
					vect[poni]<-as.numeric(le[which(le[,1]==selected[g]),2])
				}
				remain<-which(is.na(vect))
				if(length(remain)>0){
					w1<-which(interc[,1] %in% selected)
					w2<-which(interc[,2] %in% selected)
					all<-c(w1,w2)
					if(length(all)>0){
						interc<-interc[-all,]
						names2<-names2[-which(names2 %in% selected)]
						link_counts2 <-c(interc[,1],interc[,2])
						links2<-aggregate(rep(1,length(link_counts2))~link_counts2,FUN="sum")
					}
				}else{break}
			}
		}else{break}
	}
	if(length(remain)>0){
		labrem<-labels[remain]
		labcoun<-labrem
		disrem<-dis[which(rownames(dis)%in%labrem),]
		disrem[disrem == 0] <- 1000
		for(c in 1:length(remain)){
			#c<-1
			min = (which(disrem == min(disrem), arr.ind = TRUE))
			casi<-NULL
			if(!(is.null(nrow(min)))){
				casi<-nrow(min)
			}
			if(!(is.null(casi))){
				if(casi>1){
					min<-min[1,]
					na<-rownames(disrem)[min[1]]
					vect[remain[which(labrem==na)]]<-disrem[min[1],min[2]]
					disrem<-disrem[-min[1],-which(colnames(dis)==na)]
					labcoun<-labcoun[-which(labcoun==na)]
			}
				if(casi==1){
					na<-rownames(disrem)[min[1]]
					vect[remain[which(labrem==na)]]<-disrem[min[1],min[2]]
					disrem<-disrem[-min[1],-which(colnames(dis)==na)]
					labcoun<-labcoun[-which(labcoun==na)]

			}
			}
			if(is.null(casi) & length(labcoun==1)){
				vect[remain[which(labrem==labcoun)]]<-disrem[min[1]]
				
			}
			
		}
	}
	return(as.numeric(vect))
}
