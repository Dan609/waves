library(tseriesEntropy)
library(fractaldim)
library(TSclust)
library(devtools)
library(MSMVSampEn)

z <- ts(calcium$ROI422, start = calcium$Time[1])

# Entropy tests of serial and cross dependence for time series
## WARNING: computationally intensive, increase B with caution!

Srho.test.ts.p(z, B = 2, 
             lag.max = 10,
             ci.type = c("perm"),
             quant = 0.95, plot = TRUE)


Srho.ts(z, lag.max = 10, method = "integral", plot = TRUE,
        maxpts = 0, tol = 0.001)


w <- as.integer(calcium$ROI222)

ent <- Srho(w, lag.max = 10, stationary = TRUE,
     version = "FORTRAN")
ent


MSMVSampEn(mat = z, 
           M = 5, tau = 10, r = 1, eps = 1567, scaleMat = T)

print(which.max(spx))

print(which.max(spy))

print(which.max(x.spec))

print(x.spec)
