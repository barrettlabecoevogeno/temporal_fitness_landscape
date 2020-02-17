library(mgcv)
pp.tab = NULL 
ss.tab = NULL 
nn.tab = NULL 
# pdf("~/Dropbox/finch_recap/v2/data/creating_data/GCV.scores.pdf")
par(mfrow=c(4,2))
year.plot  = 2004:2010
year.plot2  = 2005:2011

for (i in 1:7) {
  # yearly.number.of.id is in the script "find local min mac.R"
mod1=model.list[[i]]
mydata=data.frame(y = yearly.number.of.id[[i]]$surv, 
                  x= yearly.number.of.id[[i]]$pc1)
mod1$gcv.ubre
my.summ=summary(mod1)
# str(my.summ)
pp.tab = c(pp.tab, list(my.summ$p.table))
ss.tab = c(ss.tab, list(my.summ$s.table))
nn.tab = c(nn.tab, list(my.summ$n))
# acf(resid(mod1))
# pacf(resid(mod1))
# plot(mod1)

sum(mod1$hat)

lambda <- exp( seq(-20,10, by=.8))        # fit a range of lambdas >0

# The GCV score is the minimised generalised cross-validation (GCV) score of the GAM fitted. 
# GCV is used for smoothness selection in the mgcv package for R; 
# smoothing parameters are chosen to minimise prediction error where ϕ is unknown, and standard CV or GCV can be used to estimate prediction error
gcvscore <- sapply(lambda, function(lambda, mydata){
  gam(y ~ s(x), data = mydata, 
      family = binomial(link = "logit"), 
      sp = lambda, method="GCV.Cp")$gcv.ubre},
  mydata)
# plot(lambda, gcvscore, type = "l");abline(h =z$gcv.ubre, lty = 3) # or

plot(log(lambda), gcvscore, type = "l", 
     main = paste(year.plot[i], 
                  year.plot2[i],sep = "-"),
     ylab ="GCV score",
     xlab = "ln(lambda)")
abline(h = mod1$gcv.ubre, v = log(mod1$sp), lty = 3)
exp.lambda = exp(mean(c(-4,-5,-4,-4,-4,-4,-4,-4,-3,-6,-13,0)))
abline(v= log(exp.lambda),col="red", lty =2)
print(nrow(mydata))
}
# dev.off()
p.tab= do.call(rbind, pp.tab)
s.tab= do.call(rbind, ss.tab)
n.tab= do.call(rbind, nn.tab)
write.csv(p.tab,"~/Dropbox/finch_recap/v2/data/creating_data/GAM.table.pp.csv")
write.csv(s.tab,"~/Dropbox/finch_recap/v2/data/creating_data/GAM.table.ss.csv")
write.csv(n.tab,"~/Dropbox/finch_recap/v2/data/creating_data/GAM.table.nn.csv")
