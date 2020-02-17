prep.data <- function(sp.keep,
                      yr.keep,
                      site.keep,
                      flip.pc1 = FALSE, # If the PCA axis is not in the same direction as before, you can flip it here 
                      adults.only = FALSE, # If true, it'll keep only females and males (remove juveniles) 
                      keep.last.year.data = TRUE, # If you don't want to drop the individuals from the last year 
                      gam.analysis = FALSE,
                      recalculate.pca = FALSE, # Recalculating the PCA for the subsetted data
                      data = NULL
                      ) {
  # load('output/bird.data.RData', verbose=TRUE)
  # source('scripts/0.1_misc.R')
  bird.data = data
  if(adults.only){
    bird.data = bird.data[bird.data$age %in% c("f","m"),]
  }
  ## reduce to relevant species
  bird.data = droplevels(bird.data[bird.data$Species1 %in% sp.keep,])
  bird.data = droplevels(bird.data[bird.data$Site %in% site.keep,])
  
  if (recalculate.pca) {
    library(vegan)
    pca.recalc = vegan::rda(bird.data[,c("MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth")],scale = FALSE)
    pc.scores = scores(pca.recalc,choices = c(1,2), scaling = 2)
    bird.data$PC1 = pc.scores$sites[,1]
    bird.data$PC2 = pc.scores$sites[,2]
  } else {pca.recalc = NULL}
  
  # If the PCA axis is not in the same direction as before, you can flip it here 
  if(flip.pc1){
    bird.data$PC1 = -1*bird.data$PC1
  }
  
  ch <- colnames(bird.data)[colnames(bird.data) %in% paste('y', yr.keep, sep='.')]
  ## remove birds that were not observed in selected years
  detected <- apply(bird.data[,ch], 1, sum)>0 # The sum in row has to be greater than 0 
  bird.data <- bird.data[detected,]
  
  ## remove birds whose first capture date was in last year (drop a lot of captures in one year)
  if(keep.last.year.data){
    bird.data <- bird.data[bird.data$first<length(ch),]
  }
  if(gam.analysis){
    bird.data[,ch] = known.state.cjs(bird.data[,ch])
    # This will only work in pairs of years (can't use it if you are doing 2003:2018, only 2003:2004, or any other combination)
    bird.data = bird.data[bird.data[,ch[1]]==1,]
  }
  ## *** individual-level variables ***
  ind.vars <- data.frame(sp      = bird.data$Species1,
                         mbl     = standardize(bird.data$MedianBeakLength),
                         mbw     = standardize(bird.data$MedianBeakWidth),
                         mbd     = standardize(bird.data$MedianBeakDepth),
                         pc1     = bird.data$PC1,
                         pc2     = bird.data$PC2,
                         band    = bird.data$BANDFINAL)
  std.mean.sd = data.frame(
    avg.scl.mbl = mean(bird.data$MedianBeakLength),
    avg.scl.mbw = mean(bird.data$MedianBeakWidth),
    avg.scl.mbd = mean(bird.data$MedianBeakDepth),
    sd.scl.mbl  = sd(bird.data$MedianBeakLength),
    sd.scl.mbw  = sd(bird.data$MedianBeakWidth),
    sd.scl.mbd  = sd(bird.data$MedianBeakDepth))
  
  list(X = as.matrix(bird.data$X.corr[,ch]), 
       first = bird.data$first, # vector of first occasion
       bir.d  = bird.data,
       ind.vars = ind.vars,
       pca.recalc = pca.recalc, 
       std.mean.sd = std.mean.sd)
}

