source("scripts/1.0_subsetting_dataset_bird.R")
plot(table(datata$Date), type = "l", lwd = .25, xaxt = "n", ylab = "Number of catches at each day")
at = seq(1,length(table(datata$Date)),by = 25)
axis(1, at = at, labels = names(table(datata$Date))[at],las =3, cex.axis = .9)

sorted.catches = sort(table(datata$Date), decreasing = TRUE)
sorted.catches.mob = sort(table(droplevels(datata[datata$By. %in% "MOB","Date"])), decreasing = TRUE)

plot(sorted.catches, 
     # xlim = c(0,100),
     type = "l", lwd = 1, xaxt = "n", ylab = "Number of catches at each day")
at = seq(1,length(sorted.catches),by = 1)
axis(1, at = at, labels = names(sorted.catches)[at],las =3, cex.axis = .8)
find.days <- function(year,names) {
  selection <- which(grepl(pattern = year,x = names(names)))
  return(selection)
}
abline(v = find.days('2018',sorted.catches), col = "red")
abline(v = find.days('2019',sorted.catches), col = "blue")
abline(v = which (names(sorted.catches) %in% names(sorted.catches.mob)), col = "yellow")
abline(h = 20, lty = 3)
