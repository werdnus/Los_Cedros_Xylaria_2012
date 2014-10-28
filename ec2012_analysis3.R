###################### community comp matrix #########################

## systematically construct a df that matches all xylaria specimens and cultures 
## to their respective name/number, grid, OTU, and life-phase (endo or strom) 
## using Roo's latest otu list, SequencesR.csv, and mushroom master spreadsheet.
## This is going to be complicated.

name_etc <- function() {
  
  library(reshape)
  setwd('/home/daniel/R/Ec2012')
  
  ##start with roo's otu sheet:
  aa <- read.csv("roo_otus_6-30-14_r.csv", as.is = TRUE)[,-4]
  names(aa)[1] <- 'OTU'
  aa <- aa[aa$cf.ID.genus == "Xylaria",]
  
  ##empty df
  bb <- as.data.frame(matrix(nrow = nrow(aa), ncol = 4))
  a <- 0
  
  ##now go cell by cell, keep all non-NAs
  for (j in 1:nrow(aa)) {
    for (i in 4:ncol(aa)) {
      
      if (is.na(aa[j,i]) == FALSE) {
        a <- a+1 ##counter
        bb[a,1] <- aa[j,i]
        bb[a,2] <- aa$cf.ID.genus[j]
        bb[a,3] <- aa$cf.ID.species[j]
        bb[a,4] <- aa$OTU[j]
      }}}
  
  colnames(bb) <- c("Name","Genus","species","OTU")
  
  ##now to assign Grids/Flag# to each Xylaria specimen/culture
  
  cc <- read.csv('SequencesR.csv', as.is = TRUE)[,1:3] ##read the spreadsheet with grid #'s, leave out seqs
  
  ##note there is an error in the SequencesR.csv, gives the grid info for 875 as 12,  
  ##really this is 12new. I will remove this here:
  
  cc[cc$Name == '875',] <- NA
  
  ##now cycle through all the xylaria specimens/cultures we got from roo's latest OTUs, 
  ## and attach Grid numbers from the SequencesR.csv spreadsheet:
  
  dd <- NULL
  
  for (i in 1:nrow(bb)) {
    
    if (any(cc$Name == bb$Name[i], na.rm = TRUE)) { 
      if (is.na(cc[which(cc$Name == bb$Name[i]),2]) == FALSE) {
        dd[i] <- cc[which(cc$Name == bb$Name[i]),2]
      }} else {dd[i] <- NA}
  }
  
  bb$Grid <- dd
  
  ##mostly works, but not finding some: '78.4.1' in the SequenceR, for instance 
  ##Not a problem, because culture names have their grid in the name, hope nothing else is missing
  
  ## attach endophyte/stromata designation
  
  ##cultures
  bb$EorS <- NA
  bb[grep('\\..\\.', bb$Name),]$EorS <- 'E' #use the periods in culture names to tell them apart
  
  ##endos everything else?
  bb[is.na(bb$EorS),]$EorS <- 'S' ####doesn't work
  
  ##filling the holes in grid info:
  ##first, check endos:
  
  bb[is.na(bb$Grid) & bb$EorS == 'E',] ##78.4.1 seems to be the only one missing, plug it in
  bb[bb$Name == '78.4.1',]$Grid <-78
  
  ##the stromata are a lot messier...
  ##read in a pruned version of mushroom master list, with just the stromata
  ## found in the plot, at a flag, not including the "new" flag we declared
  ## near plot 12, can't use this for spatial analysis. Here goes...
  
  ee <- read.csv('mush_mast_abbrev.csv', as.is = TRUE) 
  names(ee) <- c('Where', 'Grid', 'Name')
  
  ##just plots from that have numerical grid info this removes some points that just say 'BRL', 
  ##but these are also lacking grid data, so foo: 1064, 1067, 1068, 1130, all of the R collections
  ##also 's.n.' specimens these points lack grid info will be removed, what are they? Also
  ##of interest is the collection without a name/number at plot 22, found on a leaf. Doubt it is
  ##a xylaria, but I will ask. For now, it stays out. 
  
  ##so now we use this ee dframe to assign grid numbers to our stromata (better than SequencesR.csv for stromata)
  
  bb[bb$EorS == 'S',]$Grid <- NA ##clear out old Grid info for stromata, to be filled in below in only the stromata we want
  
  for (i in 1:nrow(bb)) {
    if (any(ee$Name == bb$Name[i]) == TRUE) { 
      bb$Grid[i] <- ee[ee$Name == bb$Name[i],2]
    }} 
  
  bb <- na.omit(bb) ##and get rid of specimens with NAs, should be all the ones from other years or not in our sampling scheme
  
  ##for the record, here are the OTUs that are in Roo's latest OTU sheet that are thrown out by me for 
  ##our spatial analysis:
  
  #br <- bb[bb$Name %in% ff$Name == FALSE,]
  #save(br, file = 'OTUs_removed_7-22-14.RDATA')
  
  rownames(bb) <- 1:nrow(bb)
  mgo <- bb
  mgo$OTU <- as.numeric(mgo$OTU)
  mgo$Grid <- as.integer(mgo$Grid)
  
  #save(mgo, file = "mgo.rdata")
  #write.csv(mgo, file = "mgo.csv")
  return(bb)
  
  ##and that is our new reference sheet, updated as of 5/21/2014
  ## This was a pain in the ass. goddamit, when I design my next study, specimens names will be simple numbers
}

################ spco #####################
spco2 <- function() {

  setwd("/home/daniel/R/Ec2012/")
  
  ##lets write an aggregate function for reshape, that will categorize whether or not we have
  ##all stromata, all endophytes, both, or neither.
  
  spco_code <- function(value){
    
    if (length(value) == 0) {a <- 0}else{
      if (all(value == 'S')) {a <- 1}else{
        if(all(value == 'E')) {a <- 2}else{        
          a <- 3}}}
    
    return(a)  
  }
  
  setwd('/home/daniel/R/Ec2012')  
  library(reshape)
  
  load('mgo.rdata')
  spco <- cast(mgo, Grid ~ OTU, value = 'EorS', fun.aggregate = spco_code) 
  
  ##wow, that actually worked. I don't totally understand the aggregate function, but okay.
  ##To fit into the old code for randomizations and such, fill out the zero grids
  ## and get rid of grid column:
  
  ##empty df
  aa <- mat.or.vec(nr = 120, nc = ncol(spco)-1)
  aa <- as.data.frame(aa) 
  rownames(aa) = 1:120; colnames(aa) <- colnames(spco[,-1])
  
  ##now fill with the spco.2 data
  for (i in spco$Grid) {
    aa[i,] <- spco[spco$Grid == i,-1] 
  }
  spco <- aa
  #save(spco, file = 'spco.rdata')
  #write.csv(spco, file = 'spco.csv')
  return(spco)
}

################# spatial points matrix #######################
spmat <- function(){

setwd("/home/daniel/R/Ec2012/")
library(sp)

plot_name <- 1:120

yy <- c(rep(seq(from = 2.5, to = 52.5, by = 5), times = 10),seq(from = 2.5, to = 47.5, by = (5)))   #and similar for x-axis, in multiples of 10
xx <- c(rep(seq(from = 5, to = 95, by = 10), each = 11), rep(105, times = 10))

losced_mat <- cbind(xx,yy)
row.names(losced_mat) <- plot_name

losced_pts <-as.data.frame(losced_mat)
losced_sp <- SpatialPoints(losced_mat)
summary(losced_sp)
bbox(losced_sp)

write.csv(losced_pts, file = "losced_pts.csv")
save(losced_sp, file = 'losced_sp.rdata')
save(losced_pts, file = 'losced_pts.rdata')

return(losced_pts)
}

############# endophyte dataframe ##################

EndoOnly <- function(){

setwd("/home/daniel/R/Ec2012/")
load('spco.rdata')

losced_endo <- matrix(nrow = 120, ncol = length(spco[1,]))

for (j in 1:length(spco[1,])) {
  
  for (i in 1:120) {
    
    if (spco[i,j] == 2 | spco[i,j] == 3) {losced_endo[i,j] <- 1
    } else {losced_endo[i,j] <- 0}
    
  }}

losced_endo <- as.data.frame(losced_endo)
colnames(losced_endo) <- colnames(spco); rownames(losced_endo) <- 1:120

#write.csv (losced_endo, file = "endo_pa.csv")
#save(losced_endo, file = 'losced_endo.rdata')

}
############ stromata matrix ###################

StromOnly <- function(){

setwd("/home/daniel/R/Ec2012/")
load('spco.rdata')

losced_strom <- matrix(nrow = 120, ncol = length(spco[1,]))

for (j in 1:length(spco[1,])) {
  
  for (i in 1:120) {
    
    if (spco[i,j] == 1 | spco[i,j] == 3) {losced_strom[i,j] <- 1
    } else {losced_strom[i,j] <- 0}
    
  }}

losced_strom <- as.data.frame(losced_strom)

colnames(losced_strom) <- colnames(spco); rownames(losced_strom) <- 1:120

#write.csv(losced_strom, file = "strom_pa.csv")
#save(losced_strom, file = 'losced_strom.rdata')

}

############Individual species graphics ##############

##lets make some maps!

otu_maps <- function(){

setwd("/home/daniel/R/Ec2012/")
library(sp)
load('spco.rdata')
load('losced_sp.rdata')
load('mgo.rdata') 

for (j in 1:ncol(spco)) {
  
  losced_cd.j <- data.frame(spco[,j]); colnames(losced_cd.j) = "phase"  
  losced_spcd.j <- SpatialPointsDataFrame(losced_sp, losced_cd.j, match.ID = TRUE)
  
  plcols <- c(0,"red","green","purple")
  plsymb <- c(1, 23, 22, 21)
  
  setwd("/home/daniel/R/Ec2012/sp_grids/grids_jpg")
  filename.j <- paste("OTU",colnames(spco)[j],".jpg", sep = "")
  jpeg(filename.j)
  title1.j <- paste(mgo[mgo$OTU == colnames(spco)[j],2], mgo[mgo$OTU == colnames(spco)[j],3])                    
  title2.j <- paste("OTU",colnames(spco)[j], sep = " ")
  plot(losced_spcd.j, pch = plsymb[(losced_spcd.j$phase + 1)], cex = 2, bg = plcols[(losced_spcd.j$phase + 1)])
  text(x = 55, y = 65, labels = title1.j, cex = 2)
  text(x = 25, y = -10, labels = title2.j, cex = 1.5)
  legend(x = 70, y = -2.5, leg = c("Decomposer", "Endophyte", "Both"), pch = c(23,22,21), pt.cex = 1., pt.bg = c("red","green","purple"))
  dev.off()
  
}

setwd("/home/daniel/R/Ec2012/")

}
##surprisingly, that worked.

#############species accumulation curves, adapted from Roo ##############

acc_curves <- function() {

  setwd("/home/daniel/R/Ec2012/")
  library(reshape)
  library(vegan)
  load('losced_endo.rdata')
  load('losced_strom.rdata')
  
  ##already have purely stromata and endophyte comm matrices, of just Xylaria:
  
  endoacc <- specaccum(losced_endo)
  stromacc <- specaccum(losced_strom)
  
  png('ec2012_xylaria_acc.png')
  #pdf('ec2012_xylaria_acc.pdf')
  par("mar"=c(5,5.5,3,2)+1) 
  plot(endoacc, ci.type="polygon", ci.col="#9a9a9a39", lwd=2, ylim=c(0,41), ylab="", xlab="", main="")
  par(new=T)
  plot(stromacc, ci.type="polygon", ci.col="#858585DB", lwd=2, ylim=c(0,41), ylab="", xlab="", main="") 
  title(ylab="OTUs Recovered", xlab="Plots", cex.lab=2)
  legend(x = -2, y = 42, legend=c("Endophytic Isolates","Saprotrophic Specimens"), pch=15, cex=.8, col=c("#9a9a9a39","#858585DB"))
  dev.off()
  
  specpool(losced_strom)
  ##   Species     chao  chao.se    jack1 jack1.se    jack2     boot  boot.se   n
  ##All      36 52.33333 11.70549 49.88333 4.211212 57.79958 42.18927 2.302961 120
  specpool(losced_endo)
  ##    Species chao chao.se jack1 jack1.se  jack2     boot   boot.se   n
  ##All       5    5       0 7.975 1.717617 10.925 6.105077 0.8337007 120
}

##okay, that wasn't too bad, thanks to Vegan
##now need to redo Roo's bar graph

#######################Frequency graph(s)######################

##Roo constructed a frequency plot, instead of an absolute abundance
##graph. 

##get frequency of discovery of each OTU, as endo and strom
freq_graph <- function(){
  
  setwd('/home/daniel/R/Ec2012')
  
  load('losced_strom.rdata')
  load('losced_endo.rdata')
  load('mgo.rdata')
  
  a <- colSums(losced_strom)/120
  b <- colSums(losced_endo)/120  
  
  cc <- cbind(a,b) ##make a matrix, keep their OTU#s as a names
  cc <- t(cc) ##transpose
  cc <- cc[,order(-cc[1,])] ##order it by stromata freq 
  
  ##works alright. now the names...this will be messy
  
  dd <- dimnames(cc)[[2]] ##lift the (OTU)names off of the matrix 
  ee <- NULL; ff <- 0
  gg <- as.character(mgo$OTU) ##cuz OTUs are numbers, causing funky decimals
  
  for (i in dd){
    ff <- ff+1 ##counter, since i is wonky OTU names
    ##dig which genus/species belongs to each OTU, glue them together
    ee[ff] <- paste(mgo[gg == i,][1,2], mgo[gg == i,][1,3])
  }
  
  ##so now we have a list of names that should match the order 
  ##as the OTUs appear on the graph
  
  png(file = 'specacc2.png',width = 7, height = 3.5, units = 'in', res = 700)
  #pdf(file = 'specacc2.pdf',width = 7, height = 3.5)
  mp <- barplot(cc, beside =TRUE, col=c("lightgrey","black"), names.arg = rep("", times = ncol(cc)*2),
                axes = FALSE, cex.lab = 1.25, space = c(0,.5), ylab = "Frequency")
  
  axis(2, at = c(0,.05,.1,.15,.2,.25), labels = c(0,.05,.1,.15,.2,.25), las = 2, xpd = NA)
  
  text(colMeans(mp), srt = 45, labels = ee, y = -.008, cex = .75, xpd = NA, adj = 1)
  
  legend("topright", c("Endophyte","Saprotroph"), cex=.8, 
         bty="n", fill=c("black","lightgrey"));
  dev.off()
}


########### Nearest neighbors, take 5238 or whatever #################

##for each clade, stromata points and stromphyte points will have to be considered independently, then combined 
##the flow: nearest neighbor analysis (NNA) of stromata --> NNA of stromphytes --> NNA of (stromata = stromphytes) ---- 
## --> stromphyte NNA of stromata  --> Stromata NNA of stromphyte

##with these last two, we can address the FA 
##code is messy, complicated, but I'll clean it up if we need it in the future. 


fungal_cluster <- function(nenes, cycles) {
  
  #########################Stromata NNA ######################
  
  setwd("/home/daniel/R/Ec2012/")
  library(vegan); library(sp)
  load('mgo.rdata')
  load('spco.rdata')
  load('losced_strom.rdata')
  load('losced_endo.rdata')
  load('losced_sp.rdata')
  
  pval.S_S <- NULL
  
  strom_plots <- vector("list", length(spco[1,])) #individual OTU strom spatial plots will be stored here
  
  for (i in 1:length(spco[1,])) {
    
    ## for each OTU, make a matrix of all points where that OTU is present
    strom.i <- as.matrix(losced_strom[losced_strom[,i] == 1,i]); rownames(strom.i) <- rownames(losced_strom[losced_strom[,i] == 1,])
    
    bb <- as.numeric(rownames(strom.i)) ##get the rows as numbers so as to...
    
    pts.strom.i <- losced_pts[bb,] #get spatial info for all the points with our fungus 
    
    strom_plots[[i]] <- bb #store this for downstream analyses of both endos and stroms
    
    if (sum(strom.i[,1]) < 2) {pval.S_S[i] <- NA} else {     #if there are one or less fungi in the plot in this clade, stop here
      
      dist.strom.i <- as.matrix(dist(pts.strom.i)) ##make a distance matrix of this OTU's stromata 
      
      ##now sort this, create a matrix for containing all flag points with stromata of this OTU, and the distance to their 
      ##nearest neighbors (= nenes)
      
      cc <- NULL; dd <- matrix(nrow = nrow(strom.i), ncol = nenes) 
      
      for (j in 1:length(strom.i)) {
        
        cc <- sort(dist.strom.i[,j])[2:(nenes+1)] 
        
        dd[j,] <- cc
        
      }
      
      ##for each class of neighbor (1rst, 2nd, 3rd, etc), create averages for the 
      ## entire sample area (the BRL plot) of this OTU
      
      nn_means.strom <- colMeans(dd) 
      
      ##sweet, so this seems to work for any OTU, with any number of neighbors. 
      ## at this point, OTU #i has a matrix showing the nenes of all points with stromata,
      ## an average distance in all neighbor classes (our test statistics), and it has
      ##added its list of points with stromata to a list to be used later known as 'strom_plots[[i]]'
      
      #####Stromata randomization#####
      
      nn_means.ran <- matrix(ncol = nenes, nrow = cycles) ##empty bin for 
      
      for (k in 1:cycles) {  
        
        if (length(strom.i) != 0) {     ##if there is more than one stromata for this OTU, then...
          
          random_pts <- sample(1:120, length(strom.i), replace = FALSE) ##make a random, hypethical sampling area, same # of points w/stromata
          
          random_pts_mat <- losced_pts[random_pts,] ## give these random stromata points geo info
          
          dist.ran <- as.matrix(dist(random_pts_mat)) #a distance matrix of random points
          
          ##as above, sort and 
          
          cc.ran <- NULL; dd.ran <- matrix(nrow = length(strom.i), ncol = nenes)
          
          for (j in 1:length(strom.i)) {
            
            cc.ran <- sort(dist.ran[,j])[2:(nenes+1)] 
            
            dd.ran[j,] <- cc.ran
            
          } 
          
          ##store each permutation ('cycle') result 
          nn_means.ran[k,] <- colMeans(dd.ran)} else {nn_means.ran[k,] <- NA}
      }
      
      ## so now we have a bunch of randomly generated equivalents to our single, real
      ## nearest neighbor classes from the section above. This is our distribution
      ## for our test statistics.
      
      ########## stromata p-values #################
      
      aa <- nn_means.ran
      
      for (h in 1:nenes)  {  
        
        aa <- aa[aa[,h] <= nn_means.strom[h], , drop = FALSE]      
        
      } ####this recursively subsets the randomized dataframe, to ¨nenes¨-number-of-nearest-neighbors
      
      pval.S_S.i <- length(aa)/length(nn_means.ran)
      pval.S_S[i] <- pval.S_S.i  #storing the pvalues
    }}
  
  ## so we have a p-val for this OTU, for stromata clustering to stromata
  ## if I were writing this all over again, I would break up the function here
  ## but as is, I don't want to disturb things too much.
  
  ##################### Endophyte NNA ######################
  
  ## nene number set above
  
  pval.E_E <- NULL
  endo_plots <- vector("list", length(spco[1,])) #endo spatial plots will be stored here
  
  for (i in 1:length(spco[1,])) {
    
    endo.i <- as.matrix(losced_endo[losced_endo[,i] == 1,i]); rownames(endo.i) <- rownames(losced_endo[losced_endo[,i] == 1,])
    
    #create a matrix of just endophytes from i clade
    
    bb <- as.numeric(rownames(endo.i)) 
    
    pts.endo.i <- losced_pts[bb,] #get spatial info for all the points with our fungus 
    
    endo_plots[[i]] <- bb #use for downstream analyses of both endos and stroms
    
    
    if (sum(endo.i[,1]) < 2) {pval.E_E[i] <- NA} else {     #if there are one or less fungi in the plot in this clade, stop here
      
      dist.endo.i <- as.matrix(dist(pts.endo.i))
      
      cc <- NULL; dd <- matrix(nrow = length(endo.i), ncol = nenes)
      
      for (j in 1:length(endo.i)) {
        
        cc <- sort(dist.endo.i[,j])[2:(nenes+1)] 
        
        dd[j,] <- cc
        
      }
      
      nn_means.endo <- colMeans(dd)
      
      ###sweet, so this seems to work for any clade, with any number of neighbors. 
      
      #####endo randomization#####
      
      nn_means.ran <- matrix(ncol = nenes, nrow = cycles)
      
      for (k in 1:cycles) {  
        
        if (length(endo.i) != 0) {
          
          random_pts <- sample(1:120, length(endo.i), replace = FALSE)
          
          random_pts_mat <- losced_pts[random_pts,]
          
          dist.ran <- as.matrix(dist(random_pts_mat)) #a distance matrix of random points
          
          cc.ran <- NULL; dd.ran <- matrix(nrow = length(endo.i), ncol = nenes)
          
          for (j in 1:length(endo.i)) {
            
            cc.ran <- sort(dist.ran[,j])[2:(nenes+1)] 
            
            dd.ran[j,] <- cc.ran
            
          } 
          
          nn_means.ran[k,] <- colMeans(dd.ran)} else {nn_means.ran[k,] <- NA} #averages first, second...neighbor for each randomly generate plot and stores it
      }
      
      ########## endo p-values #################
      
      aa <- nn_means.ran
      
      for (h in 1:nenes)  {  
        
        aa <- aa[aa[,h] <= nn_means.endo[h], , drop = FALSE]
        
      } ####this recursively subsets the randomized dataframe, to ¨nenes¨-number-of-nearest-neighbors
      
      pval.E_E.i <- length(aa)/length(nn_means.ran)
      pval.E_E[i] <- pval.E_E.i  #storing the pvalues
    }}
  
  ## Now I need a way to overlay the two life-history phases and analyze them from both "directions"
  ## start by overlaying the two spatial dataframes, endophyte and stromata
  
  pval.S_E <- NULL; pval.E_S <- NULL
  
  for (i in 1:length(spco[1,])) {
    
    strom_c.i <- strom_plots[[i]] #pull clade i out of lists generated by previous all-endo/all-strom NNAs
    endo_c.i <- endo_plots[[i]]
    
    if (length(strom_c.i) == 0 | length(endo_c.i) == 0) {pval.S_E[i] <- NA; pval.E_S[i] <- NA} else {
      
      combo_c.i <- strom_c.i[which (strom_c.i %in% endo_c.i == TRUE)] #shows us which pts have both endo and strom, don't need this now, maybe later.
      
      pts.strom_c.i <- losced_pts[strom_c.i,]
      pts.endo_c.i <- losced_pts[endo_c.i,]
      
      ######################### Stromata, nearest Endophyte neighbor ########################
      
      dd <- matrix(nrow = length(strom_c.i), ncol = nenes)
      
      for (k in 1:length(pts.strom_c.i[,1])) {
        
        dist.S_E.k <- as.matrix(dist(rbind(pts.strom_c.i[k,], pts.endo_c.i))) #binds one stromata point to all the endos, constructs distance matrix
        
        for (j in 1:nenes) {
          
          dd[k,(j)] <- sort(dist.S_E.k[,1])[j+1]  #extracts nearest endophyte neighbors to this one stromata
          
        }}
      
      nn_means.S_E <- colMeans(dd)
      
      
      ######################### Endophyte, nearest stromata neighbor ########################
      
      dd <- matrix(nrow = length(endo_c.i), ncol = nenes)
      
      for (k in 1:length(pts.endo_c.i[,1])) {
        
        dist.E_S.k <- as.matrix(dist(rbind(pts.endo_c.i[k,], pts.strom_c.i))) #binds one endophyte point to all the strom, constructs distance matrix
        
        for (j in 1:nenes) {
          
          dd[k,(j)] <- sort(dist.E_S.k[,1])[j+1]  #extracts nearest endophyte neighbors to this one endoata
          
        }}
      
      nn_means.E_S <- colMeans(dd)
      
      ###################################### randomization S_E (& E_S?) ######################################
      
      #####Stromata randomization##
      
      ##make a big pile of random plots with same number of stromata
      
      pts.strom_c.ran <- vector("list", length = cycles) 
      
      if (length(strom_c.i) != 0) {
        for (k in 1:cycles) {  
          
          random_pts <- sample(1:120, length(strom_c.i), replace = FALSE)
          
          pts.strom_c.ran[[k]] <- losced_pts[random_pts,]
        }}
      
      #####Endophyte randomization##
      
      ##make a big pile of random plots with same number of endophytes
      
      pts.endo_c.ran <- vector("list", length = cycles) 
      
      if (length(endo_c.i) != 0) {
        for (k in 1:cycles) {  
          
          random_pts <- sample(1:120, length(endo_c.i), replace = FALSE)
          
          pts.endo_c.ran[[k]] <- losced_pts[random_pts,]
        }}
      
      ########Stromata, randomized nearest endophyte neighbor averages############
      
      nn_means.S_E.ran <- matrix(nrow = cycles, ncol = nenes)
      
      for (j in 1:cycles) {
        
        dd.ran <- matrix(nrow = length(strom_c.i), ncol = nenes)
        
        for (l in 1:length(pts.strom_c.i[,1])) {
          
          dist.S_E.ran <- as.matrix(dist(rbind(pts.strom_c.ran[[j]][l,], pts.endo_c.ran[[j]]))) #take one point from endophyte, tack it onto all the endos, make a distance matrix
          
          for (h in 1:nenes) {
            
            dd.ran[l,h] <- sort(dist.S_E.ran[,1])[h+1]  #extracts nearest endophyte neighbors to this one stromata
            
          }}
        
        nn_means.S_E.ran[j,] <- colMeans(dd.ran)}
      
      ########Endophyte, randomized nearest stromata neighbor averages############
      
      nn_means.E_S.ran <- matrix(nrow = cycles, ncol = nenes)
      
      for (j in 1:cycles) {
        
        dd.ran <- matrix(nrow = length(endo_c.i), ncol = nenes)
        
        for (l in 1:length(pts.endo_c.i[,1])) {
          
          dist.E_S.ran <- as.matrix(dist(rbind(pts.endo_c.ran[[j]][l,], pts.strom_c.ran[[j]]))) #take one point from endophyte of randomly generated plots, tack it onto all the stroms, make a distance matrix
          
          for (h in 1:nenes) {
            
            dd.ran[l,h] <- sort(dist.E_S.ran[,1])[h+1]  #extracts nearest endophyte neighbors to this one stromata
            
          }}
        
        nn_means.E_S.ran[j,] <- colMeans(dd.ran)}
      
      ###################################### p-values S_E & E_S ########################################
      
      ########## Stromata, nearest endophyte neighbor (S_E) p-values #################
      
      aa <- nn_means.S_E.ran
      
      for (h in 1:nenes)  {  
        
        aa <- aa[aa[,h] <= nn_means.S_E[h], , drop = FALSE]
        
      } ####this recursively subsets the randomized dataframe, to ¨nenes¨-number-of-nearest-neighbors
      
      pval.S_E[i] <- length(aa)/length(nn_means.S_E.ran) #store it
      
      
      ########## Endophyte, nearest endophyte neighbor (E_S) p-values #################
      
      aa <- nn_means.E_S.ran
      
      nn_means.E_S.ran
      
      for (h in 1:nenes)  {  
        
        aa <- aa[aa[,h] <= nn_means.E_S[h], , drop = FALSE]
        
      } ####this recursively subsets the randomized dataframe, to ¨nenes¨-number-of-nearest-neighbors
      
      pval.E_S[i] <- length(aa)/length(nn_means.E_S.ran) #store it
      
    }}
  
  
  ##tie OTUs to Ju's IDs, export pvalues
  
  Ju_ID <- NULL
  aa <- 0
  
  for (k in colnames(spco)) {
    aa <- aa + 1
    Ju_ID[aa] <- paste(mgo[mgo$OTU == k,][1,2], mgo[mgo$OTU == k,][1,3])  
  }
  
  dan_NNA <- as.data.frame(cbind(Ju_ID, colnames(spco), pval.S_S, pval.E_E, pval.S_E, pval.E_S))
  colnames(dan_NNA) <- c("Species","OTU","Pval.S-S","Pval.E-E","Pval.S-E","Pval.E-S")
  save(dan_NNA, file = 'all_fung_pval.rda')
  
  ##both_lf_pval <- dan_NNA[!is.na(dan_NNA$pval.S_E),] ##gives us just the fungi with both lifecycles
  
  ##save(both_lf_pval, file = 'fungclust_pval.rda')
  
  return(dan_NNA)
  
}

################### make environmental data frame ##############################

## will construct an environmental data object that will be used for the water
## and environmental analyses

enviro_df <- function(){
  
  setwd("/home/daniel/R/Ec2012")
  
  library(sp)
  library(vegan)
  
  env <- read.csv("/home/daniel/R/Ec2012/environmental/plot_env_data.csv")
  env$eastvec <- cos(env$Aspect*pi/180) + 1; env$northvec <- sin(env$Aspect*pi/180) + 1
  env <- env[,-c(1,3,6)] ##get rid of flag# (redundant) and path info (I´m not analysing)
  
  #save(env, file = 'env.rda')
  return(env)
}
################################################ water 2.0 #################################################

## so, lets see if endophytes are tending to cluster around water
## this will mimic the stromata-around-endophyte analysis

## revising this to fit the dependencies and variable names of ec2012_analysis3.R -Dan 7/28/14

## need all of these for the water analyses functions:
setwd('/home/daniel/R/Ec2012')
load('spco.rdata')
load('losced_strom.rdata')
load('losced_endo.rdata')
load('losced_pts.rdata')
load('env.rda')

stream_geo <- function() {
  
  ##need geographical info for the stream
  env$xx <- losced_pts$xx; env$yy <- losced_pts$yy
  
  ##collapse this to just geo info of pts with a stream?
  stream_pts <- env[env$Stream == 1,][,6:7]
  
  #save(stream_pts, file = 'stream_pts.rda')
  
  return(stream_pts)
}

########## finding nene averages from data #############

##make a function that generates the average distance to nearest fungal observation of a given OTU to each point on the stream

nearest_stream <- function(which_com, nenes, OTU) {  
  
  load('stream_pts.rda')
  
  nearest_OTUs <- matrix(nrow = nrow(stream_pts), ncol = nenes) 
  ##make/reset our storage bin for distances of nearest OTU observations from individual river flags
  
  ##get geo info for OTU of interest
  ## get the geo info for points that contain our OTU:
  
  ptsOTU.i <- losced_pts[which_com[,which(colnames(which_com) == OTU)] > 0,]
  
  
  ##now we want to find the nenes of this OTU to each spot on the stream. 
  ##gotta cycle through each stream flag, here is the basic kernel:
  
  stream.j <- NULL
  
  for (j in 1:nrow(stream_pts)) {
    
    stream.j <- stream_pts[j,]; rownames(stream.j) <- "stream.j" 
    #this takes j-th row out of streams and gives a new name so rbind below won't get mad
    
    s.j.OTUs <- rbind(stream.j,ptsOTU.i) #stacks the j-th line of stream onto all of the geo info of our OTU observations
    
    dist.streamOTU.i <- as.matrix(dist(s.j.OTUs)) #constructs a distance matrix from this
    
    nearest_OTUs.j <- sort(dist.streamOTU.i[1,])[2:(nenes+1)] 
    ##sorts distances from stream to all OTU observations in ascending order, to number of nearest neighbors (+1, 1rst is self-ref)
    
    ##this needs to get stored as one row in a matrix that will hold all the nenes, columns to be averaged later: 
    nearest_OTUs[j,] <- nearest_OTUs.j
    ##and process repeated to get the nenes of other stream pts
    
  }
  
  stream_nene_avg_dist <- colMeans(nearest_OTUs) 
  ##this averages the distances of nearest observations of OTU with every flag in contact with the stream
  
  return(stream_nene_avg_dist)
  
}

############## randomizations ###############

##make a function that generates a distribution of averages distances of randomly distributed observations for a given otu. 

generate_random_dist <- function(which_com, nenes, OTU, cycles) {
  
  load('stream_pts.rda')
  
  ##nenes = number of nearest neighbors, which_com is the community matrix of choice, etc.etc.
  
  ##summing up the number of plots with a given OTU:
  
  how_many_otu.i <- sum(which_com[,which(colnames(which_com) == OTU)] > 0) 
  ran.stream_nene_avg_dist <- matrix(nrow = cycles, ncol = nenes)
  
  for (k in 1:cycles) {               #from here each cycle is one random map and nenes from it, to be stored for the dist
    
    random_pts <- sample(1:nrow(losced_pts), how_many_otu.i, replace = FALSE) #create the "observations"
    
    random_pts_mat <- losced_pts[random_pts,] #geo info of our random observations
    
    ##now we want to find the nearest random observations to each spot on the stream. 
    ##gotta cycle through each stream flag:
    
    stream.ran.j <- NULL 
    ran.nearest_OTUs <- matrix(nrow = nrow(stream_pts), ncol = nenes) #bin for the nenes of each point in the stream, for this iteration
    rownames(ran.nearest_OTUs) <- rownames(stream_pts)
    
    for (j in 1:nrow(stream_pts)) {
      
      stream.ran.j <- stream_pts[j,]; rownames(stream.ran.j) <- "stream.ran.j" 
      #this takes j-th row out of streams and gives a new name so rbind below won't get mad
      
      ran.s.j.OTUs <- rbind(stream.ran.j,random_pts_mat) #stacks the j-th line of stream onto all of the geo info of our OTU observations
      
      ran.dist.streamOTU.j <- as.matrix(dist(ran.s.j.OTUs)) #constructs a distance matrix from this
      
      ran.nearest_OTUs.j <- sort(ran.dist.streamOTU.j[1,])[2:(nenes+1)] 
      ##sorts distances from stream to all OTU observations in ascending order, to number of nearest neighbors (+1, 1rst is self-ref)
      
      ##this needs to get stored as one row in a matrix that will hold all the nenes, columns to be averaged later: 
      ran.nearest_OTUs[j,] <- ran.nearest_OTUs.j
      ##and process repeated to get the nenes of other stream pts
      
    }
    
    ran.stream_nene_avg_dist[k,] <- colMeans(ran.nearest_OTUs) }
  
  
  return(ran.stream_nene_avg_dist)}

#####################generate stream pvalues##########################

generate_stream_pval <- function(which_com, nenes, OTU, cycles) {
  
  how_many_otu.i <- sum(which_com[,which(colnames(which_com) == OTU)] > 0) 
  
  if (how_many_otu.i != 0) {
    
    nearest_stream.1 <- nearest_stream(which_com, nenes, OTU)
    generate_random_dist.1 <- generate_random_dist(which_com, nenes, OTU, cycles)
    
    aa <- generate_random_dist.1
    
    for (h in 1:nenes)  {  
      
      aa <- aa[aa[,h] <= nearest_stream.1[h], , drop = FALSE]
      
    } ####this recursively subsets the randomized dataframe, to ¨nenes¨-number-of-nearest-neighbors
    
    ##each step of the process, it retains only those with nearer h-th nearest neighbors than (or =) the data from the OTU
    
    pval.stream <- length(aa)/length(generate_random_dist.1)
    #print(aa)
    #print(nearest_stream.1)
    #print(generate_random_dist.1)
  } 
  
  else {
    
    pval.stream <- NA
  }
  
  return(pval.stream)
}

######################### sensitivity of our water pval-generator #########################

##a function to easily test out our p-val generator at different #'s of nearest neighbors

test.stream.pval <- function(which_com, nenes, OTU, cycles) {
  
  dd <- NULL  
  
  for(w in 1:nenes) {
    
    if (w > sum(which_com[,which(colnames(which_com) == OTU)] > 0)) {dd[w] <- NA} else {
      
      dd[w] <- generate_stream_pval(which_com, w, OTU, cycles)
      
    }}
  return(dd)
}

##now lets see if life_history matters
## a function to compare endos, stroms, and both 

test.stream.lifehist <- function(nenes, OTU, cycles) {
  
  aa <- test.stream.pval(losced_endo, nenes, OTU, cycles)
  bb <- test.stream.pval(losced_strom, nenes, OTU, cycles)
  cc <- test.stream.pval(spco, nenes, OTU, cycles)
  dd <- cbind(aa,bb,cc)
  
  compare <- as.data.frame(dd); colnames(compare) <- c("endo","strom","both")
  
  return(compare) 
}

############ drawing the stream ##############

##need a visual of the stream

draw_water <- function(){
  
  library(sp)
  load('losced_sp.rdata')
  load('env.rda') ##get env data
  stream <- env[,3, drop = FALSE] #keep only stream data
  stream$pch <- 21; stream$pch[which(stream$Stream == 1)] <- 22 ##change symbols
  stream$color <- "00000000"; stream$color[which(stream$Stream == 1)] <- "blue" ##and color
  
  ##set up our spatial dataframe
  losced_stream_pts <- SpatialPointsDataFrame(losced_sp, stream, match.ID = TRUE)
  
  #jpeg("water_plot.jpeg")
  plot(losced_stream_pts, pch = losced_stream_pts$pch, cex = 2, bg = losced_stream_pts$color)
  title("The path of a stream \n through Los Cedros plot")
  legend("bottomright", leg = c("Water"), pch = c(22), pt.cex = 2, pt.bg = c("blue"))
  #dev.off()
  
  save(losced_stream_pts, file = 'losced_stream_pts.rda') ##added 8/10/14 
}

################################### Environmental Analysis ################################
#################### species x habitat matrix ###########################

specieshab <- function() {
  
  setwd("/home/daniel/R/Ec2012")
  
  load('mgo.rdata')
  load('losced_pts.rdata')
  load('env.rda')
  
  ##first, need to get a observation x habitat matrix, with columns for species ID and life-stage (endo/strom) 
  ##these last two will be excluded as variables from the ordination, but used to color-code points
  ##in the final graphic
  ##we don't want to repeat flag# when one point had the same fungi multiple times, this may give the
  ##impression of habitat-specialization that is not actually present 
  ##so only one endophyte and one stromata observation of each species from each flag
  
  ##order our observations by grid, then OTU, then lifestage
  ec2012_obs <- mgo[order(mgo$Grid, mgo$OTU, mgo$EorS),] 
  
  ##an empty dataframe to begin filtering out repeats of OTUs in grids
  no_repeats <- as.data.frame(matrix(nrow = nrow(mgo), ncol = ncol(mgo)))
  colnames(no_repeats) <- colnames(ec2012_obs)
  
  no_repeats[1,] <- ec2012_obs[1,] # get the dataframe started
  
  #now write a loop that looks through the ordered dataframe, sees if the OTU above it is the same, =repeated, if so it throws this out as NA:
  
  for (i in 2:nrow(no_repeats)) {
    
    if (ec2012_obs$Grid[i] == ec2012_obs$Grid[i-1] & ec2012_obs$OTU[i] == ec2012_obs$OTU[i-1] & ec2012_obs$EorS[i] == ec2012_obs$EorS[i-1]) {  ##checks the line above for identical info
      
      no_repeats[i,] <- NA
    } else { 
      
      no_repeats[i,] <- ec2012_obs[i,]}
    
  }
  
  ##seems to work. now get rid of NAs:
  no_repeats <- na.omit(no_repeats)
  
  ##attach environmental data
  
  nmes <- c(colnames(no_repeats), colnames(env), colnames(losced_pts))
  no_repeats[,7:13] <- NA #make room for env info and geo info 
  colnames(no_repeats) <- nmes
  
  ##make a loop that matches up the three dataframes
  
  for (i in 1:nrow(no_repeats)) {
    
    no_repeats[i,7:11] <- env[which(rownames(env) == no_repeats$Grid[i]),] #find the row of env data to attach, using flag/grid#
    no_repeats[i,12:13] <- losced_pts[which(rownames(losced_pts) == no_repeats$Grid[i]),]
  }
  
  ##we need a column for proximity to water: this will be the distance of an observation 
  ##from the closest stream point. 
  
  load('stream_pts.rda')
  
  #make a small function that finds the closest stream pt to a given flag/point on our map:
  
  str_dist <- function (pt){ 
    
    load('stream_pts.rda')
    dist_mat <- rbind(losced_pts[pt,], stream_pts) #put flag# on top of stream pts
    rownames(dist_mat) <- c("pt",rownames(stream_pts)) #fix names
    dist_mat <- as.matrix(dist(dist_mat)) #make a distance matrix (Euclid)
    nearest_stream <- sort(dist_mat[1,])[2] #find the nearest stream pt
    return(nearest_stream)
    
  }
  
  for (i in 1:nrow(no_repeats)) {no_repeats$Stream[i] <- str_dist(no_repeats$Grid[i])}
  
  ##attach host data:
  
  trees <- read.csv(file = "trees_R.csv") ##tree data for each grid
  
  ##now go through the rows of the environmental dframe we're making, find the tree host
  ##from this "trees" dframe, plug it in
  
  no_repeats$host <- NA
  
  for (i in 1 : nrow(trees)) {
      
      no_repeats[no_repeats$Grid == i,]$host[] <- as.character(trees$Genus_sp[i])

     }
  
  ##get rid of specimen names and rename the dataframe:
  Spec_Hab <- no_repeats[,-1]
  
  ##make some columns factors for the analysis to come
  Spec_Hab$OTU <- as.factor(Spec_Hab$OTU)
  Spec_Hab$Slope <- as.numeric(Spec_Hab$Slope)
  Spec_Hab$EorS <- as.factor(Spec_Hab$EorS)
  Spec_Hab$host <- as.factor(Spec_Hab$host)
  
  write.csv(Spec_Hab, file = 'Spec_Hab.csv')
  save(Spec_Hab, file = 'Spec_Hab.rda')
  return (Spec_Hab)
}

###################### NMS with 5 OTUs of interest #########################

hab_NMS_2012 <- function(){
  
  library(vegan)
  library(RColorBrewer)
  library(ellipse)
  
  load('Spec_Hab.rda')
  
  combos <- c(which(Spec_Hab$OTU == "10"),which(Spec_Hab$OTU == "63"),which(Spec_Hab$OTU == "96"),which(Spec_Hab$OTU == "43"),which(Spec_Hab$OTU == "53"))
  combos_df <- Spec_Hab[combos,] ##subset to just the five OTUs that have both life-phases
  combos_df <- droplevels(combos_df) ##get the levels to match after subsetting
  
  aa <- combos_df[,6:10] ##keep just environmental data
  
  mm <- metaMDS(aa) ##run nms, stresses seem okay, found a stable solution
  
  combos_df$mds_X <- mm$points[,1]; combos_df$mds_Y <- mm$points[,2] ##add the nmds coords to combo_df 
  
  ##get some colors
  
  pal_combo <- brewer.pal(12, "Paired")[c(2,4,6,8,10)]
  
  for (i in 1:nrow(combos_df)){
    
    combos_df$Color[i] <- pal_combo[which(unique(combos_df$OTU) == combos_df[i,]$OTU)]
    
  }
  
  #save(combos_df, file = 'combos_df.rda')
  
  ##plots
  
  ##here is a general plot, without any bells/whistles
  
  #plot(mm, type = "n") ##blank plot
  #combos_df$mds_X <- mm$points[,1]; combos_df$mds_Y <- mm$points[,2] ##get the nmds coords 
  #points(cbind(combos_df$mds_X, combos_df$mds_Y), col = combos_df$Color, pch = 16) ## print points with colors
  #title(main = "Habitat specificity for five \n xylariaceous endophytes/decomposers")
  #ordihull(mm, combos_df$OTU, draw = "lines", col = "purple")
  
  ##but with this I can't color the lines according to OTU, so...
  
  ##make a function that lets us plot each OTU, one at a time
  ## will need to first plot the null plot (mm, type "n"), but 
  ## this is done in the function below, "draw_NMS"
  
  OTU_nms_plot <- function(OTU){  
    par(new = TRUE)
    aa <- combos_df[which(combos_df$OTU == OTU),]
    
    points(cbind(aa$mds_X, aa$mds_Y), col = aa$Color, pch = 16)
    #ordiellipse(m2, combos_df$OTU, show.groups = OTU, col = "green")
    ordihull(mm, combos_df$OTU, show.groups = OTU, draw = "lines", col = combos_df[which(combos_df$OTU == OTU),]$Color[1])
    
    centroid.i <- colMeans(cbind(aa$mds_X,aa$mds_Y))
    points(centroid.i[1], centroid.i[2], col = aa[which(aa$OTU == OTU),]$Color, pch = 3, cex = 2.5)
    
    
    return(aa)  
    
  }
  
  ## we can then plot the colored lines. The stats that print here are 
  ## from the perMANOVA done below.
  
  draw_NMS <- function(){
    
    png(file = "hab_spec_v3.png", width = 7, height = 7, units = 'in', res = 700)
    #pdf(file = "hab_spec_v3.pdf")  
    plot(mm, type = "n", ) 
    title(main = "Habitat specificity for five \n xylariaceous endophytes/decomposers")
    text (x = 0.42, y = -.39, "PerMANOVA F(4,57) = 1.59, P-val = 0.10", cex = .75)
    
    for (i in 1:length(unique(combos_df$OTU))) {
      
      OTU <- unique(combos_df$OTU)[i]
      
      OTU_nms_plot(OTU) 
      
    }
    ## make a legend:
    
    ots <- unique(combos_df$OTU)
    cols <- unique(combos_df$Color)
    spc <- paste('X.',unique(combos_df$species))
    
    legend("topright", legend = spc, fill = cols, cex = .8)
    
    dev.off()
  }
  
  ## this puts the perMANOVA results too far to the right on my X11/screen, 
  ## but about right on the pdf. ?
  draw_NMS()
  return(mm)
}

##################### NMS with 5 OTUs of interest by life stage #########################

## adapt above function to graph NMS of be specific to either lifestyle. Argument "EnSt" 
## is either "strom" or "endo"

## may have to adjust axis or location of labels to make the permanova results fit right,
## as is, they shift about, I don't force the scale of y or x axis here

indiv_NMS <- function(EnSt){
  
  setwd('/home/daniel/R/Ec2012/')
  library(vegan)
  library(RColorBrewer)
  library(ellipse)
  
  load('Spec_Hab.rda')
  
  combos <- c(which(Spec_Hab$OTU == "10"),which(Spec_Hab$OTU == "63"),which(Spec_Hab$OTU == "96"),which(Spec_Hab$OTU == "43"),which(Spec_Hab$OTU == "53"))
  combos_df <- Spec_Hab[combos,] ##subset to just the five OTUs that have both life-phases
  combos_df <- droplevels(combos_df) ##get the levels to match after subsetting
  
  endos_df <- combos_df[combos_df$EorS == 'E',] ##subset endo env data
  strom_df <- combos_df[combos_df$EorS == 'S',] ##subset strom env data
  
  endos_mds <- metaMDS(endos_df[,6:10]) ##run nms for endos, stresses seem okay
  strom_mds <- metaMDS(strom_df[,6:10]) ##run nms for strom, stresses seem okay
  
  endos_df$mds_X <- endos_mds$points[,1]; endos_df$mds_Y <- endos_mds$points[,2] ##add the nmds coords to df's
  strom_df$mds_X <- strom_mds$points[,1]; strom_df$mds_Y <- strom_mds$points[,2]
  
  ## assign ordination and data frame - endo or strom?
  if (EnSt == 'endo') {phase <- endos_df
                       ord <- endos_mds 
                       org <- "endophytes"                        
                       labl <- expression("PerMANOVA f(4,29)=0.45,"~R^2~"=0.06, P=0.94")
                       lx <- -0.15; ly <- -.44
  } else {
    
    if (EnSt == 'strom') {phase <- strom_df
                          ord <- strom_mds
                          org <- 'decomposers' 
                          labl <- expression("PerMANOVA f(4,24)=1.84,"~R^2~"=0.23, P=0.07")
                          lx <- -0.23; ly <- -.59   
    } else {stop("weird input: Enter either \"endo\" or \"strom\"")}}
  
  ##get some colors
  
  cols <- brewer.pal(12, "Paired")[c(2,4,6,8,10)]
  
  ##assign them
  
  for (i in 1:nrow(phase)){
    
    phase$Color[i] <- cols[which(unique(phase$OTU) == phase[i,]$OTU)]
    
  }
  
  ##make a function that lets us plot each OTU, one at a time
  ## will need to first plot the null plot (mm, type "n"), but 
  ## this is done in the function below, "draw_NMS"
  
  
  OTU_nms_plot <- function(OTU){  
    
    par(new = TRUE)
    aa <- phase[phase$OTU == OTU,]
    
    points(cbind(aa$mds_X, aa$mds_Y), col = aa$Color, pch = 16)
    ordihull(ord, phase$OTU, show.groups = OTU, draw = "lines", col = phase[phase$OTU == OTU,]$Color[1])
    centroid.i <- colMeans(cbind(aa$mds_X,aa$mds_Y))
    points(centroid.i[1], centroid.i[2], col = aa[which(aa$OTU == OTU),]$Color, pch = 3, cex = 2.0, lwd = 2.5)
    
    
    return(aa)  
    
  }
  
  ## we can then plot the colored lines. The stats that print here are 
  ## from the perMANOVA done below.
  
  draw_NMS <- function(){
    
    #png(file = paste(org,"_hab_nms.pdf", sep = ""), width = 7, height = 7, units = 'in', res = 700)
    #pdf(file = paste(org,"_hab_nms.pdf", sep = ""), width = 7, height = 7)  
    plot(ord, type = "n", ) 
    title(main = paste("Habitat specificity for five xylariaceous", org))
    text (x = lx, y = ly, labels = labl, cex = .75)
    
    for (i in 1:length(unique(combos_df$OTU))) {
      
      OTU <- unique(combos_df$OTU)[i]
      
      OTU_nms_plot(OTU) 
      
    }
    ## make a legend:
    
    ots <- unique(phase$OTU)
    cols <- unique(phase$Color)
    spc <- paste('X.',unique(phase$species))
    
    legend("bottomright", legend = spc, fill = cols, cex = .8)
    
    #dev.off()
    ## this puts the perMANOVA results too far to the right on my X11/screen, 
    ## but about right on the pdf. ?
    
  }
  
  draw_NMS()
  return(ord)
}

######################## PerMANOVA tests #############################

##use PMANOVA/NPMANOVA to test for difference among the five species, using
##the above dissimilarity matrices

perTest <- function() {
  
  setwd('/home/daniel/R/Ec2012')
  library(vegan)
  load('combos_df.rda')
  env_dat <- combos_df[,6:10] 
  
  dist_combo <- vegdist(env_dat) #BC distance of our environmental data
  
  anova(betadisper(dist_combo,combos_df$OTU))
  
  ##as I understand, betadispers calculates a sort of multi-variate variance for our groups of interest
  ##then we test with ANOVA to see if one or more variances are vastly unequal (which would be shown by a low P-val)
  ##we don't see a very low pval here (p = .30), I think it is okay to go forward
  
  a1 <- adonis(env_dat~combos_df$OTU) ##test grouping by OTU, p = .12
  a2 <- adonis(env_dat~combos_df$Stream) ##testing by distance to stream, p = .001 
  a3 <- adonis(env_dat~combos_df$EorS) ## p = .08, suggestive that endophytes and strom like diff habitat
  a4 <- adonis(env_dat~combos_df$EorS, Strata = combos_df$OTU) ##nest LifeStage in OTU, p = .08, same
  a5 <- adonis(env_dat~combos_df$OTU, Strata = combos_df$EorS) ##nest OTU in lifestage, p = same
  a6 <- adonis(env_dat~combos_df$Canopy) 
  a7 <- adonis(env_dat~combos_df$Slope)
  a8 <- adonis(env_dat~combos_df$eastvec)
  a9 <- adonis(env_dat~combos_df$northvec)

  ## lets try some linear models, multiple explanatory variables
  
  attach(combos_df)
  a10 <- adonis(env_dat~Slope+Canopy+eastvec+northvec+Stream)
  a11 <- adonis(env_dat~eastvec+Canopy+Stream+northvec+Slope) 
  a12 <- adonis(env_dat~eastvec+Slope+northvec+Stream+Canopy) 
  a13 <- adonis(env_dat~Canopy+Slope+Stream+northvec+eastvec) 
  a14 <- adonis(env_dat~Stream+Canopy+Slope+eastvec+northvec) 
  a15 <- adonis(env_dat~Stream+Canopy+Slope)
  a16 <- adonis(env_dat~Stream+Slope+Canopy)
  a17 <- adonis(env_dat~eastvec+northvec) 
  a18 <- adonis(env_dat~eastvec+Slope+northvec) 
  detach (combos_df)
  ## now lets see if there is a difference between endophytes and stromata
  
  endo1 <- combos_df[combos_df$EorS == 'E',] 
  endo2 <- env_dat[combos_df$EorS == 'E',] 
  
  strom1 <- combos_df[combos_df$EorS == 'S',] 
  strom2 <- env_dat[combos_df$EorS == 'S',] 
  
  #############look for effect of OTU in just endos/ just stromata !!!!!!!!
  #############look for clustering of all stromata, not just the 5 here, also linear models
  
  ####### Just Endophytes ############
  
  attach (endo1)
  a20 <- adonis(endo2~Slope+Canopy+eastvec+northvec+Stream) ##stream for endos: R2 = .21, P = .001
  a21 <- adonis(endo2~Canopy+eastvec+northvec+Stream+Slope) ##slope for endos: R2 = .20, p =.001
  a22 <- adonis(endo2~Slope+eastvec+northvec+Stream+Canopy) ##canopy for endos: R2 = .03, p =.015
  a23 <- adonis(endo2~Canopy+Slope+Stream+eastvec+northvec) ##aspect. little effect
  a24 <- adonis(endo2~OTU) ##OTU/endos R2=.06, P = .938
  
  ## lets check for host preference here, too. Got to rid of the the NA's,
  ## add a factor called "unknown":
  
  levels(endo1$host) <- c(levels(endo1$host), "unknown")
  endo1$host[c(2,33)] <- 'unknown'  
  a25 <- adonis(endo2~host) ##huh. very high R2 (=.619), but p = .15, probably due to vastly uneven sampling of hosts
  detach (endo1)
  
  ########### Just Stromata ##############
  
  attach (strom1)
  a26 <- adonis(strom2~Slope+Canopy+eastvec+northvec+Stream) ##stream for stroms: F=112.42, R2=.44, P=.001
  a27 <- adonis(strom2~Canopy+eastvec+northvec+Stream+Slope) ##slope for stroms: F=31.36, R2=.12, P=.001
  a28 <- adonis(strom2~Slope+eastvec+northvec+Stream+Canopy) ##canopy for stroms: F=20.61, R2=.08, p =.001
  a29 <- adonis(strom2~Canopy+Slope+Stream+eastvec+northvec) ##aspect, north. little effect 
  a30 <- adonis(strom2~Canopy+Slope+Stream+northvec+eastvec) ##aspect, east.
  a31 <- adonis(strom2~OTU) ##R2 = .24, P = .061, much more pronounced than endophytes 
  detach (strom1)
  
  
  ## doing all these tests is probably a gross abuse of perMANOVA, but I'm not sure how else to 
  ## tease out the effects of each environmental factor. This seems to be the only way to get
  ## a type III SS, by running the model with the various terms of interest at the end of the 
  ## model. I'm told this is the most conservative statement of R2 values with perMANOVA

  ## Almost all of our environmental data is significant according to this test. Compare this to 
  ## the mantel tests, which found no significant explanatory environmental variables. but I'm 
  ## confusing myself. What I am doing here is not checking to see if an environmental 
  ## variable has a direct effect on the presence or absence of species. I am merely subsetting the 
  ## environmental data on each life stage, endophyte or stromata, and ranking what factors are 
  ## most responsible for differences among the points. Low p-values are expected. Why? I am 
  ## examining the effect of multiple variables on a distance matrix constructed entirely on those 
  ## variables. Very little uncertainty in the system here.
  
  ## my understanding of adonis (PerMANOVA) is that in unbalanced designs order matters, that 
  ## collinear variables will "hog" R2 values, with variables listed first receiving more of the R2
  ## to get the most conservative estimate of variation explained in our dissimilarity 
  ## between observations that is explained by a variable, we put it after other correlated
  ## variables, or constrain permutations within "strata" arguments, (the equivalent of
  ## nesting). 
  
  ## in general, proximity to water seems to be the strongest driver of dissimilarity among 
  ## our observations, followed by slope, then canopy. Though highly significant, easterly 
  ## and northerly exposure explain less of the variation (R2 = .02), and seem to be 
  ## correlated with slope. checking correlations:
  
  cor(env_dat) 
  
  ## correlations. Nothing seems too horribly correlated. Odd that azimouth and
  ## streams are somewhat correlated, the streams at out site must run fairly
  ## straight north/south or east/west
  
  ## lets try looking at our tolerances/variance-inflation-factors (VIF). We'll 
  ## use a dummy linear model/independent variable, because vif() won't accept adonis 
  ## class objects. Doesn't matter, we're interested in the relationships among dependent 
  ## variables
  
  library(car)
    
  cc <- lm(combos_df$xx ~ Slope + Canopy + Stream + eastvec + northvec)
  vif(cc)
  
  ## all values between 1.13 and 1.5, seems okay. I guess this means my sampling is unbalanced
  ## but this not really our fault, that the enivronment in our plot is unequally different.
  ## Not much I can do, but use the most conservative SS, which I believe in adonis is done by 
  ## rerunning the model multiple times, putting term of interest last. As above. These are the 
  ## values I will report. 
  
  perPval <- list(a10,a11,a12)
  return(perPval)
}

########################Host preference questions, take II ##################################

########################## faramea endophyte preference ##########################


## Spec_Hab is a df where multiple observations of an endophyte or stromatal OTU 
## have been reduced to P/A, yes/no. Here we are going to ask if within all
## grids where our most common host, Faramea, is the host, were more of the 
## endophytes a particular endophyte than we would expect? 

## We expect that if there is no host preference, the number of times an endophyte
## appears in Faramea will be the same ratio as the number of times that endophyte
## appears in all host-trees as a group

faramea <- function(){
  
  setwd('/home/daniel/R/Ec2012')
  trees <- read.csv('trees_R.csv')
  load('Spec_Hab.rda')
  
  endos <- Spec_Hab[Spec_Hab$EorS == 'E',]
  
  ##df of endos with Faramea as host. multiple observations from one grid/tree
  ##have already been reduced to pres/abs:
  
  faramea <- na.omit(endos[endos$host ==  'Faramea_jasmin',]) 
  
  length(unique(faramea$OTU))  ##how many otus is faramea hosting? 3
  
  aa <- unique(endos$OTU) ##OTU #'s of all xylaria endophytes  
  droplevels(aa)
  
  bb <- unique(faramea$OTU) #what species of endophyte did faramea host?
  droplevels(bb)
  
  ##what proportion of the total observations in Faramea belong to each OTU?
  cc <- NULL; a <- 0
  for (i in aa) {a <- a+1; cc[a] <- sum(faramea$OTU == i)}
  rm(a)
  names(cc) <- aa
  
  ##what are the proportions of the total observations in all hosts belonging to 
  ##each OTU?
  
  dd <- NULL; a <- 0
  
  for (i in aa) {a <- a+1; dd[a] <- sum(endos$OTU == i)}
  rm(a)
  
  names(dd) <- aa
  
  ## so now we have an expected (dd), and an observed (cc)
  ## lets do a chi-squared test, using Monte Carlo simulation
  ## since we so few observations, violates assumptions of 
  ## traditional test
  
  ee <- chisq.test(cc, p = dd/sum(dd), simulate.p.value = TRUE, B = 10000)
  
  ## X-squared = 2.4492, df = NA, p-value = 0.7563
  
  ##we see no evidence that Faramea has an "endophyte preference"
  ##seems that these species are roughly present in the same ratios
  ##as we see in the plot as a whole
  
  ff <- list(cc,dd,ee)
  
  return (ff)
}

######################## X. adscendens host pref ################################

## now to test for host preference from the perspective of the endophyte, 
## using X. adscendens, our most common endophyte

## Given no host preference, Xylaria adscendens should be observed in a given 
## host with the same frequency at which that host is present  e.g. if there 
## are 3 faramea trees hosting xylariaceae in the plot and 1 danea, the faramea 
## should host three times more observations than the danea

############# all xylariaceous endophytes df ############

## for this, I will use the larger data-set of all xylariaceous endophytes,
## since our df containing just Xylaria has only five endophytes, need to
## make a new one:

allend <- function(){
  
  setwd('/home/daniel/R/Ec2012')
  trees <- read.csv('trees_R.csv')
  hOTU <- read.csv('xylariaprojectsequenses_endos_R.csv', as.is =TRUE) ##host tree, sequence Name of endos
  roo <- read.csv('roo_otus_6-30-14_r.csv', as.is = TRUE)
  
  
  ## search roo's OTU's sheet for each endo culture name & OTU#
  
  for (i in 1:nrow(hOTU)) {
    hOTU$Genus[i] <- roo[which(roo == hOTU$Name[i], arr.ind = TRUE)[1],2]
    hOTU$species[i] <- roo[which(roo == hOTU$Name[i], arr.ind = TRUE)[1],3]
    hOTU$OTU[i] <- roo[which(roo == hOTU$Name[i], arr.ind = TRUE)[1],1]
  }
  
  ##get rid of repeats
  
  ## order df:
  
  aa <- hOTU[order(hOTU$Grid, hOTU$OTU),] 
  rownames(aa) <- 1:nrow(aa)
  
  bb <- as.data.frame(matrix(nrow = nrow(aa), ncol = ncol(aa)))
  
  colnames(bb) <- colnames(aa)
  
  bb[1,] <- aa[1,] # get the dataframe started
  
  for (i in 2:nrow(aa)) {
    
    if (aa$Grid[i] == aa$Grid[i-1] & aa$OTU[i] == aa$OTU[i-1]) {  ##checks the line above for identical info
      
      bb[i,] <- NA
    } else { 
      
      bb[i,] <- aa[i,]}
    
  }
  
  ##seems to work. now get rid of NAs:
  cc <- na.omit(bb)
  
  ##remove names, meaningless now:
  dd <- cc[,-4]
  
  ##give the df a better name, save it (maybe useful later)
  all_endos <- dd
  #save(all_endos, file = 'all_endos.rda')
  return(all_endos)
}

###################### x. adsc chi-square #####################

xadsc <- function () {
  
  load('all_endos.rda') ##all endos, host and grid info
  adsc <- all_endos[all_endos$OTU == 63,] ##df of just x. adscendens
  
  ## okay, we have a df that gives pres/abs data at each grid for all
  ## endophytes, not just Xylaria.
  ## Now what?
  
  ## how many trees species are hosting endophytes? what are their 
  ## ratios?
  
  length(unique(all_endos$OTU)) ##24 spp
  
  aa <- unique(all_endos$HostID) ##hosts of all xylaria endophytes  
  
  bb <- unique(adsc$HostID) #hosts of X. adscendens
  
  ##what proportion of the total observations in X. adsc belong to each host?
  cc <- NULL; a <- 0
  for (i in aa) {a <- a+1; cc[a] <- sum(adsc$HostID == i)}
  rm(a)
  names(cc) <- aa
  
  ## what are the proportions of trees hosting any endophyte in the enire 
  ## sampling area?
  
  dd <- NULL; a <- 0
  for (i in aa) {a <- a+1; dd[a] <- sum(all_endos$HostID == i)}
  rm(a)
  names(dd) <- aa
  
  ## are the ratios of hosts of X. adsc different from the ratios of host trees in general?
  
  ee <- chisq.test(cc, p = dd/sum(dd), simulate.p.value = TRUE, B = 10000)
  
  ## X-squared = 19.7986, df = NA, p-value = 0.8623
  
  ff <- list(cc,dd,ee)
  
  return (ff)
  
}

######################## combination stream and fungal maps ############################

## now to generate maps of the 5 OTUs that displayed both life histories, over stream 
## points. These will be the most relevant maps for publication

withstream <- function() {

##lots of dependencies:

setwd("/home/daniel/R/Ec2012/")
library(sp)
load('spco.rdata')
load('losced_sp.rdata')
load('mgo.rdata') 
load('env.rda')
setwd("/home/daniel/R/Ec2012/sp_grids/str_and_fung")

## make an function that gives a (1) border and (2) fill color values  
## (3)symbols, (4) varying border thickness to indicate both the presence 
## of stream and the life stage of the fungus, use this in the loop below:

strfung <- function(zz) {
  
  aa <- NULL
  
  ##no fungi, no water:
  if (zz[1] == 0 & zz[2] ==  0) {aa[1] <- "black"; aa[2] <- "#00000000"; aa[3] <- 1; aa[4] <- 2} 
  ##no fungi, water present:
  if (zz[1] == 0 & zz[2] ==  1) {aa[1] <- "blue"; aa[2] <- "blue"; aa[3] <- 21; aa[4] <- 1}
  ##stromata, no water:
  if (zz[1] == 1 & zz[2] ==  0) {aa[1] <- "black"; aa[2] <- "red"; aa[3] <- 23; aa[4] <- 1}
  ##stromata, water present:
  if (zz[1] == 1 & zz[2] ==  1) {aa[1] <- "blue"; aa[2] <- "red"; aa[3] <- 23; aa[4] <-  4} 
  ##endophyte, no water:
  if (zz[1] == 2 & zz[2] ==  0) {aa[1] <- "black"; aa[2] <- "green"; aa[3] <- 22; aa[4] <- 1}
  ##endophyte, water present:
  if (zz[1] == 2 & zz[2] ==  1) {aa[1] <- "blue"; aa[2] <- "green"; aa[3] <- 22; aa[4] <- 4} 
  ##stromata and endophyte, no water:
  if (zz[1] == 3 & zz[2] ==  0) {aa[1] <- "black"; aa[2] <- "purple"; aa[3] <- 21; aa[4] <- 1}
  ##stromata and endophyte, water present:
  if (zz[1] == 3 & zz[2] ==  1) {aa[1] <- "blue"; aa[2] <- "purple"; aa[3] <- 21; aa[4] <- 4} 
  
  return (aa)
}

##we'll map just OTUs with both endos and stromata:
endostrom <- c('10','43','53','63','96')

for (j in endostrom){
  
  ## merge water data and OTU.j data
  losced_cd.j <- data.frame(spco[,j],env$Stream); colnames(losced_cd.j) = c("fungi", "stream")  
  
  ## give a color, symbol, and border codes, using function from above
  
  for (i in 1:nrow(losced_cd.j)){
    
    losced_cd.j$border[i] <- strfung(losced_cd.j[i,])[1]
    losced_cd.j$fill[i] <- strfung(losced_cd.j[i,])[2]
    losced_cd.j$symbol[i] <- strfung(losced_cd.j[i,])[3]
    losced_cd.j$thick[i] <- strfung(losced_cd.j[i,])[4]
  }
  
  ##convert to numbers so plot() can read them:
  losced_cd.j$symbol <- as.numeric(losced_cd.j$symbol)
  losced_cd.j$thick <- as.numeric(losced_cd.j$thick)
  
  ##make a spatial data frame
  losced_spcd.j <- SpatialPointsDataFrame(losced_sp, losced_cd.j, match.ID = TRUE)
  
  filename.j <- paste("Str_OTU",j,".pdf", sep = "")
  
  ##plot
  
  #pdf(filename.j)
  png(file = )
  plot(losced_spcd.j, pch = losced_spcd.j$symbol, cex = 2, lwd = losced_spcd.j$thick, bg = losced_spcd.j$fill, col = losced_spcd.j$border)
  
  title1.j <- paste(mgo[mgo$OTU == j,2], mgo[mgo$OTU == j,3])[1]
  title2.j <- paste("OTU", j, sep = " ") ##OTU subtitle
  
  text(x = 55, y = 65, labels = title1.j, cex = 2)
  #text(x = 25, y = -10, labels = title2.j, cex = 1.5)
  #legend(x = 70, y = -2.5, leg = c("Decomposer", "Endophyte", "Both"), pch = c(23,22,21), pt.cex = 1., pt.bg = c("red","green","purple"))
  dev.off()
  
}

setwd("/home/daniel/R/Ec2012/")

}

##################### Graph Roo's leaf-fall study ######################

## Fit leaf fall data to a linear model, create a figure

leaf <- function(){
  
  setwd("/home/daniel/R/Ec2012/")
  read.csv("2012_tree_leaf_fall.csv")->trees
  
  ## fit a model to tree 1
  
  tree1 <- function(){
    tree1 <- na.omit(trees[,1:5])  ## just tree1, now sort out leaf density by distance:
    dist <- rep(tree1$Distance,4)
    lvs <- c(tree1$Tree1.1, tree1$Tree1.2, tree1$Tree1.3, tree1$Tree1.4)
    dilv <- as.data.frame(rbind(dist, lvs))
    dilv <-dilv[ order(dist, lvs)]      
    
    aa <- as.numeric(dilv[1,]) ##dist 
    bb <- as.numeric(dilv[2,]) ##lvs
    cc <- log(bb+1) ##log lvs, have to make all values non-zero
    dd <- log(aa+1) ##log dist, have to make all values non-zero
    
    ## four models: untransformed (linear), logarithmic transformations of leaf counts
    ## or distance, both transformed log 
    
    lin <- lm(bb~aa) 
    lglv <- lm(cc~aa)
    lgdi <- lm(bb~dd)
    lgld <- lm(cc~dd)
    
    ## looks like a log decay function (really, a leaf decay function!)
    ## now draw a curve from model:
    
    dists <- seq(0,3.5,.01)
    llgdi <- predict.lm(lgdi, list(dd = dists))
    #plot(bb~aa)
    plot(bb ~ aa, xlim =c(0,25), ylim = c(0,35), col = 'red', ylab = '', xlab = '')
    lines(x = exp(dists), y = llgdi)
    
  }
  ##now tree 2, let's see if the same -log model will work
  
  tree2 <- function(){
    
    tree2 <- na.omit(trees[,-c(2:5)]) ## just tree 2
    dist <- rep(tree2$Distance,4) 
    lvs <- c(tree2$Tree2.1, tree2$Tree2.2, tree2$Tree2.3, tree2$Tree2.4)
    dilv <- as.data.frame(rbind(dist, lvs))
    dilv <-dilv[ order(dist, lvs)]      
    
    aa <- as.numeric(dilv[1,]) ##dist 
    bb <- as.numeric(dilv[2,]) ##lvs
    cc <- log(bb+1) ##log lvs, have to make all values non-zero
    dd <- log(aa+1) ##log dist, have to make all values non-zero
    
    lgdi <- lm(bb~dd) ##model, -log as with tree1 
    
    dists <- seq(0,3.5,.01)
    llgdi <- predict.lm(lgdi, list(dd = dists))
    #plot(bb~aa)  
    
    plot(bb ~ aa, xlim =c(0,25), ylim = c(0,35), col = 'blue', 
         main = 'Density of fallen leaves from two trees', pch = 5, 
         xlab = 'Distance from tree (m)', ylab =
           'leaves per meter square')
    
    lines(x = exp(dists), y = llgdi)
    
    legend(x = 8, y= 34, legend = c(expression(italic('Caryodaphnopsis theobromifolia')), 
                                    expression(italic('Garcinia madruno'))), cex=.8, bty="n", col =c("red","blue"), pch = c(1,5))
    
  }
  
  pdf('leaf_curves.pdf', width = 5, height = 5)
  aa <- tree1()
  par(new = TRUE)
  bb <- tree2()
  dev.off()
  
  cc <- list(aa,bb)
  return(cc)
}