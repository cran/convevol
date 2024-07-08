#'Plots calcConv or convSig output as a two-dimensional time series.
#'
#'plot2D Plots calcConv or convSig output as a two-dimensional time series.
#'
#'@param Ct object containing calcConvCt or convSigCt output
#'@param phy The time calibrated phylogeny of interest in phylo format
#'@param tip vector of two tip labels indicating the putatively convergent taxa to be plotted in morphospace
#'@param foc.trt vector of two traits to be used in plotting. If pca == TRUE this should indicate which PC axes should be plotted in format "PC1"
#'@param trait the matrix of trait values used to compute Ct
#'@param pca logical value indicating whether to conduct principal component analysis and plot PC scores instead of raw trait values
#'@param save logical value indicating whether an animated time series (along with individual time slices) should be saved to a user specified folder
#'@param filename character indicating the desired prefix for filenames to be saved
#'@param dir optional character indicating the folderpath of the desired save location
#'@param leg logical value indicating whether a legend should be added to plots
#'@param leg.pos character indicating the position of the legend
#'@param width pixel width of saved png files 
#'@param height pixel height of saved png files
#'@param ... optional arguments to be passed to plot
#'
#'@details None
#'
#'@return Plots tracking putatively convergent taxa in two-dimensional morphospace through time
#'@return A table with trait values, morphospace distance, and nodeheights conincident with Dmax.t. See meas.path output from calcConvCt and convSigCt 
#'
#'@import phytools
#'@import magick
#'
#'@importFrom grDevices dev.off png
#'@importFrom graphics segments
#'
#'@export
#'
#'@references Grossnickle DM, Brightly WH, Weaver LN, Stanchak KE, Roston RA, Pevsner SK, Stayton CT, Polly PD, Law CJ. 2022. A cautionary note on quantitative measures of phenotypic convergence. in revision
#'Zelditch ML, Ye J, Mitchell JS, Swiderski DL. 2017. Rare ecomorphological convergence on a complex adaptive landscape: Body size and diet mediate evolution of jaw shape in squirrels (Sciuridae). Evolution 71: 633-649
#'Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution and two new measures for quantifying and assessing the significance of convergence. Evolution 69(8): 2140-2153.
#'Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.
#'Felsenstein, J. 1985. Phylogenies and the comparative method. American Naturalist, 125, 1-15.
#'
#'@examples
#'\dontrun{
#'library(phytools)
#'library(geiger)
#'
#'# create time calibrated tree
#'mytree<-rtree(100)
#'mycalibration <- makeChronosCalib(mytree, node="root", age.max=50)
#'phy <- chronos(mytree, lambda = 1, model = "correlated", 
#'calibration = mycalibration, control = chronos.control() )
#'class(phy)<-"phylo"
#'
#'# create three normally distributed phenotypic traits
#'traits <- cbind(rnorm(Ntip(phy)),rnorm(Ntip(phy)),rnorm(Ntip(phy)))
#'colnames(traits) <- c("V1","V2","V3")
#'rownames(traits) <- phy$tip.label
#'
#'#	select two random tips, excluding sister taxa
#'pairs <- apply(combn(phy$tip.label,2),2,function(x) nodepath(phy,
#'which(phy$tip.label == x[1]),which(phy$tip.label == x[2])))
#'nosis <- combn(phy$tip.label,2)[,unlist(lapply(pairs, function(x) length(x) > 3))]
#'focaltaxa <- nosis[,sample(1:ncol(nosis),1)]
#'
#'system.time(run <- calcConvCt(phy, traits, focaltaxa))
#'system.time(run2 <- convSigCt(phy, traits, focaltaxa, nsim=100))
#'
#'plot2D(run, phy, focaltaxa[1:2], colnames(traits)[1:2], traits)
#'}

plot2D <- function(Ct, phy, tip, foc.trt, trait, pca = FALSE, save = FALSE, filename = "frame", dir = NULL, leg = FALSE, leg.pos = "topleft", width = 480, height = 480, ...){
slct <- lapply(Ct$meas.path,function(x) all(which(phy$tip %in% tip) %in% x$node))
slct.MP <- Ct$meas.path[unlist(slct)][[1]]
slct.MP <- slct.MP[order(slct.MP$height),]

trt.indx <- which(colnames(trait) %in% foc.trt)
trt.name <- colnames(slct.MP)[grepl("node.anc",colnames(slct.MP))][trt.indx]

if(pca == TRUE){
PC <- prcomp(trait)

node.tmp <- slct.MP[,grepl("node.",colnames(slct.MP))]
colnames(node.tmp) <- colnames(trait)
slct.MP[,grepl("node.",colnames(slct.MP))] <- predict(PC,node.tmp)

int.tmp <- slct.MP[,grepl("int.",colnames(slct.MP))]
colnames(int.tmp) <- colnames(trait)
slct.MP[,grepl("int.",colnames(slct.MP))] <- predict(PC,int.tmp)


trt.indx <- which(colnames(PC$x) %in% foc.trt)
trt.name <- colnames(slct.MP)[grepl("node.anc",colnames(slct.MP))][trt.indx]
}

Dmx <- which(slct.MP$anc.diff == max(na.omit(slct.MP[!slct.MP$path == "tip",]$anc.diff)))
if(!is.null(Ct[["limits"]])){
	slct.cons <- slct.MP[slct.MP$height < Ct$limits[which(slct == TRUE)],]
	Dmx <- which(slct.cons$anc.diff == max(na.omit(slct.cons[!slct.cons$path == "tip",]$anc.diff)))
	}


if(!is.null(dir)){setwd(dir)}

for(n in 1:(nrow(slct.MP)-1)){
foc.df <- slct.MP[1:n,]

P1.plot <- 
apply(foc.df[!foc.df$path == "tip",],1,function(x){if(x["path"] %in% c("P1","mrca")){as.numeric(x[trt.name])} else {as.numeric(x[gsub("node","int",trt.name)])}
							})

P2.plot <- 
apply(foc.df[!foc.df$path == "tip",],1,function(x){if(x["path"] %in% c("P2","mrca")){as.numeric(x[trt.name])} else {as.numeric(x[gsub("node","int",trt.name)])}
							})

if(n >= nrow(slct.MP)-1){
	tpz <- slct.MP[slct.MP$path == "tip",]
	tpz <- tpz[order(as.numeric(tpz$node)),]
	P1.plot <- cbind(P1.plot,"tip" = as.numeric(tpz[1,trt.name]))
	P2.plot <- cbind(P2.plot,"tip" = as.numeric(tpz[2,trt.name]))
}

	if(save == TRUE){png(width = width, height = height, paste(paste(filename,10+n,sep = "_"),".png",sep = ""))}
if(pca == TRUE) {
plot(PC$x[,foc.trt], col = "grey80", asp = 1, ...)
} else {
plot(trait[,foc.trt], col = "grey80", asp = 1, ...)
}
points(t(P1.plot)[,1],t(P1.plot)[,2], type = "l", col = "grey")
points(t(P1.plot)[,1],t(P1.plot)[,2], pch = 16, col = c(rep("grey",n-1),"black"), cex = c(rep(0.75,n-1),1))

points(t(P2.plot)[,1],t(P2.plot)[,2], type = "l", pch = 16, col = "grey")
points(t(P2.plot)[,1],t(P2.plot)[,2], pch = 16, col = c(rep("grey",n-1),"black"), cex = c(rep(0.75,n-1),1))

if(n >= Dmx) {
DMX <- rbind(t(P1.plot)[Dmx,],t(P2.plot)[Dmx,])
points(DMX, pch = 8)
segments(DMX[1,1],DMX[1,2],DMX[2,1],DMX[2,2],lty = 3)
}

if(leg == TRUE){
legend(leg.pos,pch = c(1,16,8),lty = c(0,1,3),legend = c("observed tips","focal lineage","Dmax.t"),col = c("grey80","black","black"),merge = FALSE)
}
	if(save == TRUE){dev.off()}
}

if(save == TRUE){
gif.img <- lapply(list.files(pattern = filename),image_read)
gif <- image_animate(image_join(gif.img),fps = 2)
image_write(gif,paste(filename,"animation.gif",sep = "_"))
}

if(save == TRUE){
m <- nrow(slct.MP)-1

if(pca == TRUE) {
plot(PC$x[,foc.trt], col = "grey80", asp = 1, ...)
} else {
plot(trait[,foc.trt], col = "grey80", asp = 1, ...)
}
points(t(P1.plot)[,1],t(P1.plot)[,2], type = "l", col = "grey")
points(t(P1.plot)[,1],t(P1.plot)[,2], pch = 16, col = c(rep("grey",m-1),"black"), cex = c(rep(0.75,m-1),1))

points(t(P2.plot)[,1],t(P2.plot)[,2], type = "l", pch = 16, col = "grey")
points(t(P2.plot)[,1],t(P2.plot)[,2], pch = 16, col = c(rep("grey",m-1),"black"), cex = c(rep(0.75,m-1),1))

DMX <- rbind(t(P1.plot)[Dmx,],t(P2.plot)[Dmx,])
points(DMX, pch = 8)
segments(DMX[1,1],DMX[1,2],DMX[2,1],DMX[2,2],lty = 3)
if(leg == TRUE){
legend(leg.pos,pch = c(1,16,8),lty = c(0,1,3),legend = c("observed tips","focal lineage","Dmax.t"),col = c("grey80","black","black"),merge = FALSE)
}
}

return(slct.MP[Dmx,])
}
