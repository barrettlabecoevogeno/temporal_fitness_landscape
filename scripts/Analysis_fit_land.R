# # # # # # # # # # # # # # # # # # # # 
# Author: Marc-Olivier Beausoleil 
# Date: January 17, 2019
# McGill University 
# Analysis of the dynamics of Darwin's finches' fitness landscapes  
# # # # # # # # # # # # # # # # # # # # 


# Preparation of variables and data  --------------------------------------
source('scripts/0.0_initialize.R')
load('output/bird.data.RData', verbose=TRUE)
table(bird.data$Species1,
      bird.data$Year)

# Select the variables to run the scripts ---------------------------------
save.data = "./output/biotic.factors.on.survival_Andrew_meeting_changing_PCA_SCORE_for_only_fortis.RData"

sp.list <- c("fortis", "fuliginosa", "scandens","magnirostris")
yr = yr.list.subset[6:length(yr.list.subset)] # yr.list.subset[c(7,8,9,10,11,12,13,15,16)]
jump = 1
site.list <- "El Garrapatero"

pdf.output <- FALSE
find.peaks.and.valleys = FALSE
standard.ized = TRUE
# If you want only the linear βX not the βx and γX^2
linear.only = TRUE  
# If you want orthogonal X values: BEWARE! This is not the same as modeling an ORTHOGONAL regression 
orthogonal.x = FALSE

# prepare the variables that I'm going to record while iterating 
fit.grad.table = list()
model.list.logistic = NULL
pseudr = NULL
gofitfit = NULL
glm.model.list = NULL
yearly.number.of.id = NULL
all.ranges.x = NULL
list.min.fit = NULL
list.min.fit.trait = NULL
original.x = NULL
midd.list = NULL
my.eco.evo.df = NULL
full.data = NULL
final.df.new.analysis = NULL
newlist = list()
oldxlist = list()
oldzlist = list()
old_beak_L_list = list()
old_beak_W_list = list()
old_beak_D_list = list()
old_band_list = list()
model.list = list()
survived.list = list()

# These were the relatively good values of lambda used 
exp.lambda = exp(mean(c(-4)))

if(pdf.output){
  pdf(paste0("output/my.fit.land.short2_jump",
             jump,"_",gsub('([[:punct:]])|\\s+','_',site.list),".pdf"),
      height = 7,
      width = 8)
}

with(bird.data, table(Species1,Year))
# Function finding maximum and minimum of the fitness function  -----------
par(mfrow=c(1,1))
# For loop that will calculate the GAM, show the landscape and 
# will let you select what is the maximum and minimum of the function 
plot.gam.cust <- function(mod,bss=NULL,kk,title ="", bdr = TRUE) {
  kk.t=ifelse(kk==-1,"default",kk)
  vis.gam(mod,view=c("X1","X2"), main ="",
          color="heat",n.grid=50, type="link", plot.type="contour", nCol=50)
  title(main = paste(title,"Type:",bss,"&",kk.t,"dimensions,","sp=",smooth.par,sep = " "),
        cex.main = .8, col.main= "black")
  points(x = mydata$X1, y = mydata$X2, pch = 21, bg = mydata$sp, col = mydata$sp)
  points(x = mydata[mydata$y %in% 1 ,"X1"], y = mydata[mydata$y %in% 1 ,"X2"], pch = 21, bg = "yellow", col = "yellow", cex = .7) # plot only the one that survived
  vis.gam(mod,view=c("X1","X2"), main = "3D view",
          theta=40,phi=40,color="heat",n.grid=50, ticktype="detailed",type="link", plot.type="persp", border = bdr)
}

plot.gcv.score1 <- function(mod,data,bss,kk) {
  lambda <- exp( seq(-20,10, by=.8))        # fit a range of lambdas >0
  gcvscore <- sapply(lambda, function(lambda, data){
    gam(y~s(X1, bs = bss, k = kk)+s(X2, bs = bss, k = kk), 
        family = binomial(link = "logit"),
        sp = lambda,
        data=data, method="GCV.Cp")$gcv.ubre},
    mydata)
  
  plot(log(lambda), gcvscore, type = "l", 
       main = paste(yr.list[1], 
                    yr.list[2],sep = "-"),
       ylab ="GCV score",
       xlab = "ln(lambda)")
  abline(h = mod$gcv.ubre, 
         v = log(mod$full.sp[1]),
         lty = 3)
}

plot.gcv.score2 <- function(mod,data, bss,kk) {
  lambda <- exp( seq(-20,10, by=.8))        # fit a range of lambdas >0
  
  gcvscore <- sapply(lambda, function(lambda, data){
    gam(y~s(X1, bs = bss, k = kk) + s(X2, bs = bss, k = kk) + s(X1,X2, bs = bss, k = kk.i), 
        family = binomial(link = "logit"),
        sp = lambda,
        data=mydata, method="GCV.Cp")$gcv.ubre},
    mydata)
  
  plot(log(lambda), gcvscore, type = "l", 
       main = paste(yr.list[1], 
                    yr.list[2],sep = "-"),
       ylab ="GCV score",
       xlab = "ln(lambda)")
  abline(h = mod$gcv.ubre, 
         v = log(mod$full.sp[1]),
         lty = 3)
}

smooth.par = -7
kk = 10
# bss = "tp"  #c("tp","ts","ds","cr","cc","ps","cp", ### maybe
# "cs", "sos","re","mrf","gp","so","sw","sf") ### NOPEEEE 
# ts?
# ds?
# tp?
# CR?
# CC?
# CP?
bssss = c("tp")#,"ts","ds") # "ps","cp","cc","cr", ()
# if(find.peaks.and.valleys){
raw.data = NULL
model.list1 = NULL
model.list2 = NULL
plot(bird.data$PC2~bird.data$PC1, pch =".")
text(bird.data[bird.data$PC2 >.5,"PC1"], bird.data[bird.data$PC2 >.5,"PC2"], 
     labels = bird.data[bird.data$PC2 >.5,"BANDFINAL"], cex = .5)

pdf("~/Desktop/my.fit.test.pdf")
for (j in 1:length(bssss)) {
  bss = bssss[j]

  for(i in 1:c(length(yr)-1)){
    par(mfrow=c(3,2))
    
    yr.list <- c(yr[i],yr[i+jump])
    
    mdat <- prep.data(sp.keep = sp.list,
                      yr.keep = yr.list,
                      site.keep = site.list,
                      flip.pc1 = FALSE,
                      adults.only = FALSE, # If true, it'll keep only females and males (remove juveniles) 
                      keep.last.year.data = FALSE,
                      gam.analysis = TRUE,
                      recalculate.pca = FALSE,
                      data = bird.data) 

    # This is to calcualte the GAM 
    # getting response variable and the explanatory variable 
    y = as.vector(mdat$X[,2])
    x = cbind(mdat$ind.vars$pc1,mdat$ind.vars$pc2) # Works for everything execpt 2007-2008
    mbd = c(mdat$ind.vars$mbd) 
    mbl = c(mdat$ind.vars$mbl) 
    mbw = c(mdat$ind.vars$mbw) 
    band = as.character(mdat$ind.vars$band) 
    year.var = rep(yr.list[2],length(mdat$ind.vars$pc1))
    sp.dat = mdat$ind.vars$sp
    mydata = data.frame(x,y,sp = sp.dat)
    dat.for.comparison.analysis = data.frame(x,y,year.var)
    full.data = c(full.data,list(dat.for.comparison.analysis))
    # plot(mydata$X2,mydata$X1)
    if (length(which(table(mydata$sp,mydata$y)[,2] == 0)) > 0) {
      rm.sp.name = names(which(table(mydata$sp,mydata$y)[,2] == 0))
      mydata = droplevels(mydata[!(mydata$sp %in% rm.sp.name),])
    }
if (length(which(mydata$y == 1)) <5) {
  next
}

    raw.data = c(raw.data,list(mydata))
    categorical_interact1 = NULL
    categorical_interact2 = NULL
    categorical_interact3 = NULL
    
    
    cat("\n ---------- Iteration",i,"----------\n\n")
    
    # Model 1 
    categorical_interact1 <- gam(y~s(X1, bs = bss, k = kk) + s(X2, bs = bss, k = kk),
                                sp = exp(rep(smooth.par,2)),
                                data=mydata, 
                                family = binomial(link = "logit"))
    
    # Model 2 
    nb.par = length(unique(mydata$sp))
    categorical_interact1.1 <- gam(y~s(X1, bs = bss, k = kk, by = sp) + s(X2, bs = bss, k = kk, by = sp),
                                sp = exp(rep(smooth.par,2*nb.par)),
                                data=mydata, 
                                family = binomial(link = "logit"))
    
    kk.i = ifelse(kk==-1,-1,kk+2)
    # Model 3 
    categorical_interact2 <- gam(y~ s(X1, X2, bs = bss, k = kk, by = sp),
                                sp = exp(rep(smooth.par,2+nb.par)),
                                data=mydata, 
                                family = binomial(link = "logit"))
    model.list1 = c(model.list1, list(categorical_interact1))
    model.list2 = c(model.list2, list(categorical_interact2))
    
    print(summary(categorical_interact1))
    print(summary(categorical_interact1.1))
    print(summary(categorical_interact2))
    summary(categorical_interact3)
    
    # plot the models 
    plot.gam.cust(mod = categorical_interact1,bdr = TRUE, bss = bss,kk = kk,title = "GAM PC1-2,")
    plot.gam.cust(mod = categorical_interact1.1, bss = bss,kk = kk,title = "GAM by sp, PC1-2,")
    # FIND GCV SCORE 
    # plot.gcv.score1(mod = categorical_interact1,data = mydata)
    
    plot.gam.cust(mod = categorical_interact2, bss = bss,kk = kk,title = "Gam interac. PC1-2,")
    # FIND GCV SCORE 
    # plot.gcv.score2(mod = categorical_interact2,data = mydata)
    
    # plot.gam.cust(mod = categorical_interact3, bss = bss,title = "Tensor product smooths PC1 and 2")
    title(paste("Fitness landscape in year",paste(yr.list, collapse = " ")), line = -1.5, outer = TRUE)
  }
}
dev.off()

aic1 = lapply(model.list1, function(x) x$aic); which.min(unlist(aic1)); min(unlist(aic1))
aic2 = lapply(model.list2, function(x) x$aic); which.min(unlist(aic2)); min(unlist(aic2))
p.val.mod1.1 = lapply(model.list1, function(x) summary(x)$s.table[1,"p-value"])
p.val.mod1.2 = lapply(model.list1, function(x) summary(x)$s.table[2,"p-value"])
which(unlist(p.val.mod1.1) < 0.05)
which(unlist(p.val.mod1.2) < 0.05)
p.val.mod2.1 = lapply(model.list2, function(x) summary(x)$s.table[1,"p-value"])
p.val.mod2.2 = lapply(model.list2, function(x) summary(x)$s.table[2,"p-value"])
p.val.mod2.3 = lapply(model.list2, function(x) summary(x)$s.table[3,"p-value"])
which(unlist(p.val.mod2.1) < 0.05)
which(unlist(p.val.mod2.2) < 0.05)
which(unlist(p.val.mod2.3) < 0.05)
#Approximate significance of smooth terms:
#         edf Ref.df Chi.sq p-value  
#  s(X1) 6.710  7.771 17.213  0.0208 *
#  s(X2) 7.156  7.960  9.359  0.3322  

    # open3d()
    # plot3d(x = mydata$X1, y = mydata$X2, mydata$y, #type="n",
    #        xlab="PC1", ylab="PC2", zlab="Apparent Survival", 
    #        main =paste("Fitness landscape",yr.list[2], sep =" "),
    #        axes=TRUE, box=TRUE, aspect=1,col=ifelse(mydata$y, "orange", "blue"), 
    #        size = 11)
  
    # vis.gam(categorical_interact,view=c("X1","X2"),theta=40,# n.grid=500,
    #         n.grid=50,border=NA)
    # vis.gam(categorical_interact3,view=c("X1","X2"),theta=40,# n.grid=500,
    #         n.grid=50)
    # library(visreg)
    # visreg(categorical_interact)
    # library(mgcViz)
    # car::scatter3d(y~X1+X2, data=mydata,fit=c("additive"),df.additive	 = 5)
    
    # Fitting the GAM 
    z <- gam(y ~ s(X1), 
             data = mydata, 
             family = binomial(link = "logit"), # Logistic link 
             # sp = exp.lambda,
             method = "GCV.Cp")
    
    model.list = c(model.list,list(z))
    
    yr1 = substr(i + 2003, 3, 4)
    yr2 = substr(i + 2003 + jump, 3, 4)
    oldxlist = c(oldxlist, list(x))
    oldzlist = c(oldzlist, list(z$fitted.values))
    old_beak_L_list = c(old_beak_L_list, list(mbd))
    old_beak_W_list = c(old_beak_W_list, list(mbl))
    old_beak_D_list = c(old_beak_D_list, list(mbw))
    old_band_list = c(old_band_list, list(band))
    
    survived.list = c(survived.list, list(y))
    
    lambb.z = round(log(exp.lambda))
    
    # Getting new x that is spaced evenly respecting the GAM function (fitness function) 
    newx <- seq(from = min(mydata$X1), 
                to = max(mydata$X1), 
                length.out = 2000)
    
    # Using the model to generate the new response variable 
    z1 <- predict(z, 
                  newdata=list(X1 = newx), 
                  se.fit = TRUE)
    
    # This is the actual transformed data 
    yhat <- invlogit(z1$fit)
    upper <- invlogit(z1$fit + z1$se.fit)
    lower <- invlogit(z1$fit - z1$se.fit)
    
    # Here is the plot of the fitness function using the evenly spaced 
    # data (newx) and the response to it (yhat)
    plot(newx, yhat, type="l", 
         ylim = c(0,1), 
         xlab = "PC1",
         main = bquote(atop("GAM±1SE",
                            lambda * " = "*.(lambb.z) * 
                              ",  yr = "*.(yr1) *"-"* .(yr2))))
    # Adding error 
    lines(newx, upper, lty = 2)
    lines(newx, lower, lty = 2)
  }
    # find the *minimum* of the fitness function from the gam by clicking on BOTH sides of the highest visible peak and  the minimum value between the 2 peaks (valley) in the GAM  
    # midd  = locator(n = 2)
    # Starting from the left side, click right around the maximum of the fitness function. This will find the maximum value 
    # peak  = locator(n = 4)
    # midd=NULL
    # This will extract the X values for the selections 
    # midd = midd$x
    # peak = peak$x
    # midd = readline(prompt="Enter middle to find local minima: ")
    
    # Here is the actual function that will find the maximum and minium 
    # local.min  = min(yhat[newx > midd[1] & newx < midd[2]])
    # local.max  = max(yhat[newx > midd[1] & newx < midd[2]])
    
    # make a database usign the 2 peaks 
    # local.max.peak1  = newx[which(yhat==max(yhat[newx > peak[1] & newx < peak[2]]))]
    # local.max.peak2  = newx[which(yhat==max(yhat[newx > peak[3] & newx < peak[4]]))]
    
    # Draw a diagnostic line to see that you've selected the mximium of an obserbed peak and the minimum of the valley 
    # abline(h = c(local.min, local.max), lty = 3)
    # This is showing the 2 maximum values of the 2 peaks 
    # abline(v = c(local.max.peak1, local.max.peak2), lty = 3, col = "red")
    
    # This si a metric of the fitness function itself. That means that the highest fitness- the lowest fitness would be the fitness differential (if it exists)
    # midd.list = c(local.max - local.min)
    
    # Make a record of all the expected response varaible from the evenly spaced X 
    # newlist = c(newlist, list(yhat))
    
    # Make a database for all these new varaibles 
    my.eco.evo.df=rbind(my.eco.evo.df,data.frame(mid = midd.list,
                                                 # yr = paste(yr.list,collapse = '_'), 
                                                 yr1 = yr.list[1],
                                                 yr2 = yr.list[2], 
                                                 local.max.peak1 = local.max.peak1,
                                                 local.max.peak2 = local.max.peak2,
                                                 preci.yr1 = mdat$yr.vars$precipitation.year.no.std[1], 
                                                 preci.yr2 = mdat$yr.vars$precipitation.year.no.std[2],
                                                 sum.preci.yr1 = mdat$yr.vars$sum.rainfall.no.std[1], 
                                                 sum.preci.yr2 = mdat$yr.vars$sum.rainfall.no.std[2]))
    
  }# End of for(i in 1:c(length(yr)-1)){
  
  eff = structure(c(36, 140, 212, 120, 52, 56, 
                    132, 300, 128, 120, 128, 
                    104, 30.2565802161513, 81.4673085182818, 
                    81.4673085182818, 117.180771285717
  ), .Names = c("2003", "2004", "2005", 
                "2006", "2007", "2008", 
                "2009", "2010", "2011", 
                "2012", "2013", "2014", 
                "2015", "2016", "2017", "2018"))
  my.eff = eff[names(eff)%in% yr]
  
  par(mfrow=c(2,2))
  heights = my.eco.evo.df[,1]
  my.eco.evo.df$good = c(1,1,0,1,1,1,0,0,0,0,0,0,0)
  my.eco.evo.df$good = c(1,1,1,1,1,1,1,0,0,0,0,0,0)
  my.eco.evo.df$good = rep(1, nrow(my.eco.evo.df))
  my.eco.evo.df$peak.height = as.factor(c("l","r","r","r",
                                          "l","r","l","r",
                                          "l","m","m","m",
                                          "m"))
  my.eco.evo.df$two.peaks = c(0,1,1,1,
                              1,1,1,0,
                              0,0,0,0,
                              0)
  my.eco.evo.df[my.eco.evo.df$peak.height == "m",""]
  # Keep only the relevant year combination
  my.eco.evo.df$effort = my.eff[1:c(length(my.eff)-1)]
  my.eco.evo.df$effort = my.eco.evo.df$effort/max(my.eco.evo.df$effort)
  my.df = my.eco.evo.df[my.eco.evo.df$good == 1,]
  
  fct = (mid)~(sum.preci.yr2)
  lm.out1 = lm(fct,data = my.df)
  summary(lm.out1)
  plot(fct,data = my.df, 
       cex = my.df$effort, 
       pch =21, 
       bg = my.df$peak.height, 
       col = my.df$peak.height,
       ylab = "Valley depth (in fitness difference)",
       xlab = "Cumulative precipitation");abline(lm.out1)
  text(x = my.df$sum.preci.yr1,
       y = my.df$mid,
       labels = my.df$yr2, cex = .6, pos = 3)
  
  fct = log((mid))~log((sum.preci.yr2))
  lm.out2=lm(fct,data = my.df)
  summary(lm.out2)
  plot(fct, data = my.df, 
       cex = my.df$effort, 
       pch =21, 
       bg = my.df$peak.height, 
       col = my.df$peak.height, 
       ylab = "ln(Valley depth) (in fitness difference)",
       xlab = "ln(Cumulative precipitation)");abline(lm.out2)
  text(x = log(my.df$sum.preci.yr1),
       y = log(my.df$mid),
       labels = my.df$yr2, cex = .6, pos = 3)
} # End of if(find.peaks.and.valleys){


# Save finding peaks ------------------------------------------------------
if(find.peaks.and.valleys){
  save(my.eco.evo.df,
       lm.out1, lm.out2,
       newlist,
       newx,
       model.list,
       oldxlist, oldzlist,
       old_beak_L_list, old_beak_W_list, old_beak_D_list,
       old_band_list,
       survived.list,
       full.data,
       file = save.data)
}