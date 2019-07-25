# ROI01 save to z0 for coherence calibration

calcium <- read.csv2('094.csv', sep=',')
colnames(calcium) <- c("frame", 
                       'ROI01', 'ROI02', 'ROI03', 'ROI04', 'ROI05', 
                       'ROI06', 'ROI07', 'ROI08', 'ROI09', 'ROI10',
                       'ROI11', 'ROI12', 'ROI13', 'ROI14', 'ROI15', 'ROI16')
calcium$time <- calcium$frame * 110.803 / 1567
z0 <- ts(calcium$ROI01, start = calcium$time[1])
tail(calcium)

# choose variable
# 'ROI_512x512', 'ROI_180x40', 'ROI_40x40', 'ROI_10x10', 'ROI_2x2'
# Rename columns
calcium <- read.csv2('094px10.csv', sep=',')

colnames(calcium) <- c("frame", 
                       'ROI_512x512', 'ROI_180x40', 'ROI_40x40', 'ROI_10x10')


head(calcium)

# Rename columns
# calcium <- read.csv2('094px.csv', sep=',')
# colnames(calcium) <- c("frame", 
#                       'ROI_512x512', 'ROI_180x40', 'ROI_40x40', 'ROI_2x2', 'ROI_10x10', 'ROI_2x2_2')
#------------------------------------------
# Calibrate time: 110.803 sec / 1567 frames

calcium$time <- calcium$frame * 110.803 / 1567

tail(calcium)

#------------------
# sinus calibration

calcium$sin <- c(rep(NA, 1567))

fs <- 14.142211 # samples per second

df <- 1/fs # seconds per second

stoptime <-  110.803 # seconds

t <- seq(0, stoptime, by = df)

fsin <- 1

sinwave <- sin(2*pi*fsin*t)

plot(t, sinwave, type = 'l',
     main = "sinus wave",
     xlab = 'Time, sec', ylab = 'value',
     xlim = c(0, 10), ylim = c(-1, 1)
)

calcium$sin <- sinwave

xcal <- calcium$sin




# time series analysis

calcium_ts <- ts(xcal, start = calcium$time[1])

del <- 0.0707103
x.spec <- spectrum(calcium_ts, log = "no", span = 10, plot = FALSE)
spx <- x.spec$freq/del
spy <- 1*x.spec$spec

plot(spy ~ spx, xlab="frequency (Hz)", ylab="spectral density",
     type="l", lwd = 2)



# Frequency

z <- ts(calcium$sin, start = calcium$time[1])

spectrogram(z, 14.29745, windowlength = 95000, 0.01,
            timestep = 1000,
            maxfreq = 0.5,
            colors = TRUE,
            dynamicrange = 40, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'hamming', 
            quality = TRUE)

# period

xcal <- z

x <- ts(xcal, start = calcium$time[1])

my.data <- data.frame(x = x)

my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/500,
                        lowerPeriod = 2,
                        upperPeriod = 64,
                        make.pval = TRUE, n.sim = 10)

wt.image(my.w, n.levels = 250)


# wavelet coherency

x <- z
y <- z0

# y <- z

my.data <- data.frame(x = x, y = y, date = calcium$time)

my.wc <- analyze.coherency(my.data, my.pair = c("x", "y"),
                           loess.span = 0,
                           dt = 1, 
                           dj = 1/100,
                           lowerPeriod = 4, upperPeriod = 32,
                           make.pval = TRUE, n.sim = 30)


wc.image(my.wc, n.levels = 250, color.key = "interval",
         siglvl.contour = 0.1, siglvl.arrow = 0.05, which.arrow.sig = "wt",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "Time",
         label.time.axis = TRUE,
         main = "")





# peaks detection

spec(xcal, f = 14.142211, nmax = 5,
     threshold = 100,
     plot = T,
     title = TRUE,
     #xlab = "", ylab = "",
     labels = TRUE, legend = TRUE, digits = 5)












acf(xcal, lag.max = 100, main = "Autocorrelation")

# plot results

par(mfrow=c(2,3))

plot(xcal ~ calcium$time, 
     xlab = "time (s)", ylab = "fluorescence intensity (a.u.)", t="l", 
     xlim = c(0, 1), ylim = c(0, 1), main = "Total signal")

plot(xcal ~ calcium$time, 
     xlab = "time (s)", ylab = "fluorescence intensity (a.u.)", t="l", 
     xlim = c(0, 20), ylim = c(0, 1), main = "Signal segment")

boxplot(xcal, ylab = "fluorescence intensity (a.u.)",
        main = "Mean intensity", ylim = c(0, 150))

plot(spy ~ spx, xlab="frequency (Hz)", ylab="spectral density",
     type="l", lwd = 2)

acf(xcal, lag.max = 500, main = "Autocorrelation")



# calibration

z <- ts(calcium$ROI_512x512, start = calcium$time[1])

z <- ts(calcium$ROI_180x40, start = calcium$time[1])

z <- ts(calcium$ROI_40x40, start = calcium$time[1])

z <- ts(calcium$ROI_10x10, start = calcium$time[1])


fs <- 14.142211 # samples per second

df <- 1/fs # seconds per second

stoptime <-  110.803 # seconds

t <- seq(0, stoptime, by = df)

plot(z ~ t, type = 'l', main = "calcium$ROI",
     xlab = 'Time, sec', ylab = 'Ca',
     xlim = c(0, 110), ylim = c(0, 2000)
     )

fsin <- 1

sinwave <- sin(2*pi*fsin*t)

plot(t, sinwave, type = 'l',
     main = "sinus wave",
     xlab = 'Time, sec', ylab = 'value',
     xlim = c(0, 10), ylim = c(-1, 1)
)


tail(calcium)

tail(sinwave)

calcium$sin <- c(rep(NA, 1567))

calcium$sin <- sinwave
