# path to directory where package is being developed
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata/'
path3 = '/slm/data/'

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
bou[,"trt"] <- as.factor(bou[,"trt"])
bou[,"y"] <- 7-bou[,"i"]
bou[,"x"] <- bou[,"j"]
caribouDF = bou
# write object to disk as a data object for the package
path3 = 'data/'
save(caribouDF, file = paste0(path1, path3, 'caribouDF.rda'))

