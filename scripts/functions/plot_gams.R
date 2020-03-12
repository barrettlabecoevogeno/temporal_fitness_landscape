# # # # # # # # # # # # # # # # # # # # 
# Author: Marc-Olivier Beausoleil 
# Date: 12 March 2020
# McGill University 
# Functions to plot the gams 
# # # # # # # # # # # # # # # # # # # # 


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
