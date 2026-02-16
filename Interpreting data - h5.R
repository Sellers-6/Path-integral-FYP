# nolint start

library(ggplot2)
library(gridExtra)  # side by side plots
library(dplyr)
library(rhdf5)

h5ls("data.h5")

# vec <- h5read("dat3.h5", "/E0Therm/Periodic/QHO")
# head(vec)

# bc <- "Periodic"    # or "Dirichlet"
bc <- "Dirichlet"   # or "Periodic"
sys <- "QHO"        # or "DWP"
# sys <- "DWP"        # or "QHO"

# Read data

thermSweeps <- as.numeric(unlist(h5read("data.h5", paste0("/thermSweeps/", bc, "/", sys))))

E0ThermData <- as.numeric(unlist(h5read("data.h5", paste0("/E0Therm/", bc, "/", sys))))

accRateThermData <- as.numeric(unlist(h5read("data.h5", paste0("/accRateTherm/", bc, "/", sys))))

E0DecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/E0Decorr/", bc, "/", sys))))

accRateDecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/accRateDecorr/", bc, "/", sys))))

psiDecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/psiDecorr/", bc, "/", sys))))

# Variables

a <- 2              # DWP variables
lambda <- 1/12

measures <- 200

repeats <- length(E0DecorrData) / measures

pathLength <- 5000

#### Thermalisation ####

# Some quick sanity checks

min(thermSweeps)  # Not too small (<1000)
max(thermSweeps)  # Not too big (>200000)


mean(accRateThermData) * 100 # Should be ~ 60 - 80%

# Track E0 as it thermalises, to best see this, an average of E0 for the first "minThermSweeps" is found

thermalisationInterval <- 100
thermMeasures <- ceiling(thermSweeps / thermalisationInterval)
sum(thermMeasures) == length(E0ThermData) # Needs to be True

E0Groups <- rep(seq_along(thermMeasures), thermMeasures)

E0Split <- split(E0ThermData, E0Groups)

minLen <- min(sapply(E0Split, length))
E0Trimmed <- lapply(E0Split, `[`, 1:minLen)

E0Mat <- do.call(cbind, E0Trimmed)
E0Avg <- rowMeans(E0Mat)


ggplot(data.frame(index = 1:length(E0Avg), E0 = E0Avg), aes(x = index * thermalisationInterval, y = E0)) +
    geom_point() +
    labs(
      x = "Measurement index",
      y = "E0",
      title = paste("Thermalisation of", bc, sys)
    )

# Decorrelation

E0DecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/E0Decorr/", bc, "/", sys))))

E0Split <- split(E0DecorrData, rep(1:repeats, each = measures)) # Split into list of repeats

E0RepeatAvg <- sapply(E0Split, mean)  # Take mean per repeat

length(E0RepeatAvg)
bins <- 13

hE0 <- hist(E0RepeatAvg, breaks = bins, plot = FALSE) # Compute histogram breaks and bin width
binWidth <- diff(hE0$breaks)[1]

E0Norm <- hE0$counts / (sum(hE0$counts) * binWidth) # Normalised histogram counts (probability density)

continuousE0 <- seq(min(E0RepeatAvg), max(E0RepeatAvg), length.out = 200) # Continuous x values for normal curve
normDist <- dnorm(continuousE0, mean = mean(E0RepeatAvg), sd = sd(E0RepeatAvg))
dx <- diff(continuousE0)[1]
normDist <- normDist / sum(normDist * dx)  # Normalisation

# Histogram and normal curve
histPlot <- ggplot() +
  geom_histogram(aes(x = E0RepeatAvg, y = after_stat(density)),
                 bins = bins, fill = "skyblue", color = "black") +
  geom_line(aes(x = continuousE0, y = normDist),
            color = "red", linewidth = 1) +
  labs(title = "Ground State Energy Histogram",
       x = "Energy", y = "Probability Density") +
  theme_minimal()

# QQ plot
shapiroTest <- shapiro.test(E0RepeatAvg)
qqPlot <- ggplot(data.frame(E0RepeatAvg), aes(sample = E0RepeatAvg)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = paste0("QQ Plot (Shapiro-Wilk p = ", round(shapiroTest$p.value, 4), ")"),
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

# Combined plots side by side
grid.arrange(histPlot, qqPlot, ncol = 2)

histPlot

mean(E0RepeatAvg)
sd(E0RepeatAvg)

# Wave function

# matrixOfPositions <- as.numeric(as.matrix(waveFunction[[name]])) 
# xVec <- as.vector(matrixOfPositions) 
xVec <- as.numeric(unlist(h5read("data.h5", paste0("/psiDecorr/", bc, "/", sys))))

bins <- 100

hist(xVec, breaks = bins, main = "Histogram of positions",
     xlab = "position", ylab = "count")
h <- hist(xVec, breaks = 100, plot = FALSE)

counts <- sum(h$counts)                             # Total counts which equals measures * lattice size
positionRange <- range(h$mids)[2] - range(h$mids)[1]  # Range of positions 
binWidth <- positionRange / bins                      # Width of each bin 
normFactor <- counts * binWidth                     # Normalization factor 
psi <- sqrt(h$counts / normFactor)                    # Normalized wave function 

if (sys == "QHO") {
  psiAnalytical <- exp(-(h$mids ^ 2) / 2)               # The analytical wavefunction of the QHO
} else if (sys == "DWP") {
   x0 <- sqrt(1 / (2 * a * lambda))   # approximate minima positions
 
  # Double Gaussian superposition
  psiAnalytical <- exp(-((h$mids - x0)^2)/2) +
                   exp(-((h$mids + x0)^2)/2)
}

psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) # Normalising the analytical wavefunction 

ggplot(data.frame(x = h$mids, psi = psi), aes(x = x, y = psi)) +
  geom_line() +
  geom_point() +
  geom_line(aes(x = h$mids, y = psiAnalytical)) +
  labs(title = "Wave Function", x = "Position", y = "Psi")

# Correlation function

xVec <- as.numeric(unlist(h5read("data.h5", paste0("/psiDecorr/", bc, "/", sys))))

expectedLength <- repeats * measures * pathLength
length(xVec) == expectedLength # Ensure xVec is the correct length

groupSize <- measures * pathLength

repeatIDs <- rep(1:repeats, each = groupSize)

xVecSplit <- split(xVec, repeatIDs)

xVecSplitMat <- lapply(xVecSplit, function(v) {
  matrix(v, nrow = pathLength, ncol = measures)
})


twoPointCorrelator <- function(xMat) {
  N <- nrow(xMat)  # path length
  nMeas <- ncol(xMat)
  
  # Initialize correlation vector
  corr <- numeric(N)
  
  # Loop over path positions
  for (t in 1:N) {
    # For each offset n = 0:(N-1)
    for (n in 0:(N-1)) {
      # Circular indexing like C++ modulo
      idx <- ((t - 1 + n) %% N) + 1
      # Average over measurements
      corr[n + 1] <- corr[n + 1] + mean(xMat[t, ] * xMat[idx, ])
    }
  }
  
  # Normalize by path length
  corr <- corr / N
  return(corr)
}

twoPointCorrelator_fast <- function(xMat) {
  # Average over measurements first
  xAvg <- rowMeans(xMat)
  N <- length(xAvg)
  # Circular convolution via FFT
  corr <- Re(fft(Conj(fft(xAvg)) * fft(xAvg), inverse = TRUE)) / N
  return(corr)
}

correlationList <- lapply(xVecSplitMat, twoPointCorrelator_fast)

# correlationList <- lapply(xVecSplitMat, twoPointCorrelator)

correlationAvg <- Reduce("+", correlationList) / length(correlationList)

# ---- Plot ----
dfCorr <- data.frame(
  lag = 0:(length(correlationAvg)-1),
  correlation = correlationAvg
)

ggplot(dfCorr, aes(x = lag, y = correlation)) +
  geom_line(color = "#000000") +
  labs(
    title = paste("Decorrelation function for", bc, sys),
    x = "Time (index of the path)",
    y = "G(t, 0)"
  )

# E1 from the two-point correlator



# nolint end
