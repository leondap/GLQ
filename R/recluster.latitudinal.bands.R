recluster.latitudinal.bands<-function (location, identification, sequences, coordinates, populations, range = c(0, 90), width = 3, asympt = FALSE, mini_sp = 2) {
    interval <- width/2
    start <- range[1]
    slices <- round(((range[2] - range[1])/interval) + 1, 0)
    res <- NULL
    matrixsp <- array(NA, dim = c(1e+05, 17))
    dataset <- cbind(location, identification, coordinates, populations)
    species <- unique(identification)
    matricerichness <- array(0, dim = c(length(species), slices + 
        1))
    rownames(matricerichness) <- species
    matricehaplotypes <- array(0, dim = c(length(species), slices + 
        1))
    rownames(matricehaplotypes) <- species
    matricehaplotypeslow <- array(0, dim = c(length(species), 
        slices + 1))
    rownames(matricehaplotypeslow) <- species
    matricedistances <- array(NA, dim = c(length(species), slices + 
        1))
    rownames(matricedistances) <- species
    matriceGD <- array(NA, dim = c(length(species), slices + 
        1))
    rownames(matriceGD) <- species
    matriceGDred <- array(NA, dim = c(length(species), slices + 
        1))
    rownames(matriceGD) <- species
    matriceST <- array(NA, dim = c(length(species), slices + 
        1))
    rownames(matriceGD) <- species
    matriceSTred <- array(NA, dim = c(length(species), slices + 
        1))
    rownames(matriceGD) <- species
    matriceALL <- array(NA, dim = c(length(species), slices + 
        1))
    rownames(matriceGD) <- species
    matriceALLred <- array(NA, dim = c(length(species), slices + 
        1))
    rownames(matriceGD) <- species
    matricedistancesGD <- array(NA, dim = c(length(species), 
        slices + 1))
    rownames(matricedistancesGD) <- species
    matricedistancesST <- array(NA, dim = c(length(species), 
        slices + 1))
    rownames(matricedistancesST) <- species
        for (cut in 0:slices) {
#cut<-20
            lati <- start + interval * cut
            stripe <- which(dataset[, 1] >= lati - interval & 
                dataset[, 1] < lati + interval)
            if (length(stripe) > 0) {
                datalat <- dataset[stripe, ]
                speciessp <- aggregate(rep(1, nrow(datalat)) ~ datalat[, 2], FUN = "sum")                  
                for (k in 1:nrow(speciessp)) {
                  spe <- which(rownames(matricerichness) == speciessp[k, 1])                    
                  matricerichness[spe, cut + 1] <- speciessp[k,  2]
                   
                }
            }
        }
        line <- 1
        for (spe in 1:length(species)) {
#spe<-259
            dataspeq <- dataset[which(dataset[, 2] == species[spe]), ]
            seqspeq <- sequences[which(dataset[, 2] == species[spe])]
            stripes <- which(matricerichness[spe, ] >= mini_sp)
            if (length(stripes) > 0) {
                minimum <- min(matricerichness[spe, ][stripes])
                max <- max(matricerichness[spe, ])
                for (strip in 1:length(stripes)) {
#strip<-8
                  lati <- start + interval * (stripes[strip] - 1)
                  stripe <- which(dataspeq[, 1] >= lati - interval & dataspeq[, 1] < lati + interval)
                  seqslat <- seqspeq[stripe]
                  datalat <- dataspeq[stripe, ]
                  baricenter <- c(mean(datalat[, 3]), mean(datalat[, 4]))
                  coordin <- rbind(baricenter, datalat[, 3:4])
			qualiseqs <- sample(1:length(seqslat))[1:minimum]
                  molt <- 1
                  same <- as.vector(dist(datalat$populations))
			popul<-(rep(datalat$populations,length(datalat)))
                  same[which(same > 0)] <- 1
                  regi <- unique(datalat$populations)
                  regi <- regi[!is.na(regi)]
			numberregi<-length(regi)
                  GDdis <- NULL
                  numero <- 0
                  for (re in 1:length(regi)) {
                    piglia <- datalat[which(datalat$populations == regi[re]), ]                      
                    if (nrow(piglia) > 1) {
                      numero <- numero + nrow(piglia)
                      baricenterGD <- c(mean(piglia[, 3]), mean(piglia[, 4]))
                      coordinGD <- rbind(baricenterGD, piglia[, 3:4])                        
                      distanzGD <- as.matrix(earth.dist(coordinGD, dist = TRUE))[1, ]
                      distanzGD <- distanzGD[2:length(distanzGD)]
                      GDdis <- c(GDdis, distanzGD)
                    }
                  }
                  molt<-1
                  tutte <- dist.dna(seqslat, model = "raw", pairwise.deletion = TRUE)
                  tuttered <- dist.dna(seqslat[qualiseqs], model = "raw", pairwise.deletion = TRUE)
                  matriceALL[spe, stripes[strip]] <- mean(tutte) 
                  matriceALLred[spe, stripes[strip]] <- mean(tuttered)
                  intra <- which(same == 0)
                  matriceGD[spe, stripes[strip]] <- mean(as.vector(tutte)[intra]) 
                  distanz <- as.matrix(earth.dist(coordin, dist = TRUE))[1, ]
                 distanz <- distanz[2:length(distanz)]
                 matricedistancesGD[spe, stripes[strip]] <- sqrt((sum(GDdis^2)/numero))
                 inter <- which(same == 1)
                 matriceST[spe, stripes[strip]] <- mean(as.vector(tutte)[inter])
                 matricehaplotypeslow[spe, stripes[strip]] <- length(recluster.haplotypes(seqslat[qualiseqs])$frequency) 
                 matricedistances[spe, stripes[strip]] <- sqrt(mean(distanz^2))
                 matrixsp[line, 1] <- species[spe]
                 matrixsp[line, 2] <- matricedistances[spe, stripes[strip]]
                 matrixsp[line, 3] <- matriceALL[spe, stripes[strip]]
                 matrixsp[line, 4] <- matricehaplotypeslow[spe, stripes[strip]]                    
                 matrixsp[line, 11] <- lati
                 matrixsp[line, 7] <- length(seqslat)
                 matrixsp[line, 8] <- length(qualiseqs)
                 matrixsp[line, 9] <- baricenter[1]
                 matrixsp[line, 10] <- baricenter[2]
                 matrixsp[line, 6] <- matriceALLred[spe, stripes[strip]]
                 matrixsp[line, 12] <- matricedistancesGD[spe, stripes[strip]]                    
                 matrixsp[line, 13] <- matriceGD[spe, stripes[strip]]                 
                 matrixsp[line, 14] <- matriceST[spe, stripes[strip]]
			matrixsp[line, 15] <- numberregi
			matrixsp[line, 16] <- length(inter)
			matrixsp[line, 17] <- length(intra)


                  if (asympt) {
                    if (length(freq) > 1) {
                      accu <- iNEXT(freq, q = 0, datatype = "abundance", 
                        size = NULL, endpoint = max, knots = 40, 
                        se = TRUE, conf = 0.95, nboot = 50)
                      matricehaplotypes[spe, stripes[strip], 
                        it] <- accu$AsyEst[1, 2] * molt
                      matrixsp[line, 5] <- matricehaplotypes[spe, 
                        stripes[strip]]
                    }
                    if (length(freq) == 1) {
                      matricehaplotypes[spe, stripes[strip]] <- 1 * molt
                      matrixsp[line, 5] <- matricehaplotypes[spe, 
                        stripes[strip]]
                    }
                  }
                  line <- line + 1
                }
            }
        }
    colnames(matrixsp) <- c("species", "Distances", "IGV", "haplotypes_red", 
        "haplotype_asy", "IGV_red", "Specimens", "Specimens_red", 
        "Bar_long", "Bar_lat", "slice", "Distances_GD", "IGV_intra", 
        "IGV_inter", "regions", "Cells_inter","Cells_intra")
    rowgood <- which(complete.cases(matrixsp[, 1:2]))
        res$separate <- as.data.frame(matrixsp[rowgood,])
       return(res)
}
