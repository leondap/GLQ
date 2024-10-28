recluster.landscape<-function (mat, units = NULL, dist, transcorr=null, map_alt, map1=NULL, mapplot=NULL, minsea = 2, plot=F) 
{
	gendistances <- as.matrix(dist)
	mat[, 1] <- mat[, 1] + runif(nrow(mat), min = 0.001, max = 0.002)
	mat[, 2] <- mat[, 2] + runif(nrow(mat), min = 0.001, max = 0.002)
	try <- deldir(mat[, 1], mat[, 2])
	maxdat <- max(max(try$delsgs[, 5]), max(try$delsgs[, 6]))
	if (maxdat < nrow(mat)) {
		stop("Not all objects in triangulation consider collinearity")
	}
	midpointx <- NULL
	midpointy <- NULL
	gendist <- NULL
	lengthpath <- NULL
	lengthsea <- NULL
	seapointx <- NULL
	seapointy <- NULL
	altimin <- NULL
	altimax <- NULL
	if(plot){
		plot(mat)
		if(!(is.null(mapplot))){
			sp::plot(mapplot, add=T)
		}
	}
	for (i in 1:nrow(try$delsgs)) {
#i<-2
		gendist[i] <- gendistances[try$delsgs[i, 5], try$delsgs[i, 6]]
		if (abs(try$delsgs[i, 1] - try$delsgs[i, 3]) < 0.043 & abs(try$delsgs[i, 2] - try$delsgs[i, 4]) < 0.043) {
			try$delsgs[i, 1] <- try$delsgs[i, 1] + ((0.043 - abs(try$delsgs[i, 1] - try$delsgs[i, 3])) * sign(try$delsgs[i, 1] - try$delsgs[i, 3]))
			try$delsgs[i, 2] <- try$delsgs[i, 2] + ((0.043 - abs(try$delsgs[i, 2] - try$delsgs[i, 4])) * sign(try$delsgs[i, 2] - try$delsgs[i, 4]))
		}
		if(!(is.null(transcorr))){
			patch <- gdistance::shortestPath(transcorr, c(try$delsgs[i, 1], try$delsgs[i, 2]), c(try$delsgs[i, 3], try$delsgs[i, 4]), output = "SpatialLines")
			pa <- as(patch, "sf")
			st_crs(pa) <- 4326
		}else{
			punto1 <- st_point(c(try$delsgs[i, 1], try$delsgs[i, 2]))  # Punto 1: long=10.0, lat=45.0
			punto2 <- st_point(c(try$delsgs[i, 3], try$delsgs[i, 4]))  # Punto 2: long=12.0, lat=47.0
			geom1 <- st_sfc(punto1, crs = 4326)  # CRS WGS84 (EPSG:4326)
			geom2 <- st_sfc(punto2, crs = 4326)
			pa <- st_sfc(st_cast(st_union(geom1, geom2), "LINESTRING"))
		}
		if(plot){
			plot(pa, add=T)
		}
      	nume <- as.numeric(round(st_length(pa) * 0.0003, 0))
		punti <- st_sample(pa, nume, type = "regular")
		coords <- (st_coordinates(punti)[, 1:2])
		coords <- rbind(c(10, 10), coords)
		value <- c(as.vector(extract(map1, coords)))
		altitudes <- as.data.frame(as.vector(extract(map_alt, coords)))
		val <- cbind(coords, value, altitudes)
		altimin[i] <- min(val[, 4], na.rm = T)
		altimax[i] <- max(val[, 4], na.rm = T)
		lengthpath[i] <- st_length(pa)*0.001
		if (nume == 1) {
			midpointx[i] <- (try$delsgs[i, 1] + try$delsgs[i, 3])/2
            	midpointy[i] <- (try$delsgs[i, 2] + try$delsgs[i, 4])/2
            	seapointx[i] <- (try$delsgs[i, 1] + try$delsgs[i, 3])/2
            	seapointy[i] <- (try$delsgs[i, 2] + try$delsgs[i, 4])/2
            	lengthsea[i] <- 0
		}
		if (nume > 1) {
			val <- val[-1, ]
			tabe <- with(rle(val[, 3]), {
			ok <- values == 1
			ends <- cumsum(lengths)[ok]
			starts <- ends - lengths[ok] + 1
			resi <- as.data.frame(cbind(starts, ends))
			return(resi)
		})
		tabe <- cbind(tabe, (tabe[, 2] - tabe[, 1] + 1))
		tabe <- tabe[which(tabe[, 3] >= minsea), ]
		line <- which.max(tabe[, 3])[1]
		midpointx[i] <- as.numeric(val[round(nrow(val)/2, 0), 1])
		midpointy[i] <- as.numeric(val[round(nrow(val)/2, 0), 2])
		lengthsea[i] <- sum(tabe[, 3])/3
		if (lengthsea[i] > 0) {
            	seapointx[i] <- as.numeric(val[round((tabe[line, 1] + tabe[line, 2])/2, 0), 1])
			seapointy[i] <- as.numeric(val[round((tabe[line, 1] + tabe[line, 2])/2, 0), 2])
		}
		if (lengthsea[i] == 0) {
			seapointx[i] <- midpointx[i]
			seapointy[i] <- midpointy[i]
			}
		}
	}
	table <- cbind(try$delsgs, gendist, midpointx, midpointy, 
	lengthpath, seapointx, seapointy, lengthsea, altimin, altimax)
	return(table)
}
