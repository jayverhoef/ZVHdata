library(slm)
library(maptools)
library(lattice)
library(spdep)
library(nlme)
library(classInt)
library(viridis)

# ---------- CREATE THE DATA FRAME 

bou <- data.frame(
  i = as.vector(t(outer(1:6, rep(1, times = 5)))),
  j = as.vector(outer(1:5, rep(1, times = 6))),
  water = c("Y", "Y", "Y", "N", "N",
    "Y", "N", "Y", "Y", "Y",
    "Y", "Y", "Y", "N", "N",
    "N", "N", "N", "Y", "Y",
    "N", "Y", "N", "N", "N",
    "N", "N", "N", "Y", "Y"),
  tarp = c("clear", "shade", "none", "clear", "shade",
    "none", "clear", "clear", "shade", "none",
    "clear", "shade", "none", "shade", "none",
	  "clear", "shade", "none", "shade", "none",
	  "none", "clear", "clear", "shade", "none",
	  "clear", "shade", "none", "clear", "shade"),
  trt = c( 2,  1,  3,  5,  4,
           3,  5,  2,  1,  3,
           2,  1,  3,  4,  6,
           5,  4,  6,  1,  3,
           6,  2,  5,  4,  6,
           5,  4,  6,  2,  1),
  z = c(2.421, 2.443, 1.810, 1.966, 2.377,
        2.224, 2.101, 1.804, 1.963, 2.099,
        1.886, 2.130, 1.861, 2.366, 1.994,
        1.637, 2.447, 1.819, 2.045, 1.841,
        2.115, 1.996, 1.801, 2.049, 2.129,
        2.023, 2.050, 2.029, 1.789, 2.063)
)
bou[,"trt"] <- as.factor(data[,"trt"])
bou[,"y"] <- 7-data[,"i"]
bou[,"x"] <- data[,"j"]

cip = classIntervals(bou$z, n = 4, style = 'fisher')
palp = viridis(6)[3:6]
cip_colors = findColours(cip, palp)

pdf('/media/jay/data/DaleBook/Chapter1/figure/caribou_Design.pdf', 
  width = 10, height= 8)
layout(matrix(1:2, nrow = 1), widths = c(4,1))
par(mar = c(5,5,1,1))
plot(bou$x,bou$y, pch =15, cex = 12, xlim = c(0.7,5.3), ylim = c(0.7,6.3),
  col = cip_colors, cex.lab = 2.5, cex.axis = 2, xlab = 'x-coordinate',
  ylab = 'y-coordinate')
text(bou$x,bou$y,labels = bou$tarp, pos = 3, cex = 2)
text(bou$x,bou$y,labels = bou$water, pos = 1, cex = 2)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .2, xright = .3, ytop = .8,
	breaks = cip$brks, colors = palp, cex = 2.5)
dev.off()


pdf('/media/jay/data/DaleBook/Chapter1/figure/caribou_boxplots.pdf', 
  width = 10, height= 7)
par(mar = c(2,6,1,1))
boxplot(z ~ tarp + water, data = bou, cex.lab = 2.3, cex.axis = 1.5,
  ylab = 'Percent Nitrogen in Prostrate Willows', pch = 19, cex = 2, lwd = 2)
dev.off()


