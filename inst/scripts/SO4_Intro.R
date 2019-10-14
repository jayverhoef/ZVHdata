library(slm)
#library(sp)
library(classInt)
library(viridis)
library(colorspace)

# load data from package
data(SO4obs)
data(USboundary)

# Make traditional plots of predictions and standard errors
cip = classIntervals(SO4obs@data$SO4, n = 6, style = 'fisher')
palp = viridis(6)
cip_colors = findColours(cip, palp)

pdf('/media/jay/data/Zimmerman Book/Chapter1/figure/SO4.pdf', width = 11, height= 8.5)
layout(matrix(1:2, nrow = 1, byrow = TRUE), widths = c(3,1))
par(mar = c(0,0,5,0))
plot(SO4obs, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
	breaks = cip$brks, colors = palp, cex = 1.5)
dev.off()

system('pdfcrop /media/jay/data/DaleBook/Chapter1/figure/SO4.pdf
  /media/jay/data/DaleBook/Chapter1/figure/SO4-crop.pdf')
