# nolint start

library(ggplot2)
library(gridExtra)  # side by side plots
library(dplyr)
library(rhdf5)

# Variables from the simulation, needed for interpreting the data

a <- 3              # DWP variables
lambda <- 0.1

S_inst <- sqrt(lambda / 3) * (2 * (a ^ 3)) / 3 # Two equivelant ways of calculating the instanton action
S_inst
S_inst <- (2 / 3) * (a ^ 2) 
S_inst

alpha <- 1 / 12 # A complicated calculation performed in Zinn-Justin 1993 gives this value

K <- sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5) # A prefactor for the splitting energy
K
4 / sqrt(pi)

Splitting_energy <- K * exp(-S_inst)
Splitting_energy

E0_inst <- 0.5 - (Splitting_energy / 2)
E1_inst <- 0.5 + (Splitting_energy / 2)
E0_inst
E1_inst

measures <- 100

repeats <- 30

pathLength <- 5000

latticeSpacing <- 0.1
thermalisationInterval <- 10000

acceptableError <- 0.01 # This was the ratio of the monte carlo error in ground state energy to the current average ground state energy

# Boundary conditions 

# bc <- "Periodic"
bc <- "Dirichlet"   # I have noticed that Dirichlet is systematically worse than periodic boundary conditions

# System type

# sys <- "QHO"
# sys <- "AHO"
sys <- "DWP"

# Read data

thermSweeps <- as.numeric(unlist(h5read("data.h5", paste0("/thermSweeps/", bc, "/", sys))))

E0ThermData <- as.numeric(unlist(h5read("data.h5", paste0("/E0Therm/", bc, "/", sys))))

accRateThermData <- as.numeric(unlist(h5read("data.h5", paste0("/accRateTherm/", bc, "/", sys))))

xThermData <- as.numeric(unlist(h5read("data.h5", paste0("/xTherm/", bc, "/", sys))))

E0DecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/E0Decorr/", bc, "/", sys))))

accRateDecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/accRateDecorr/", bc, "/", sys))))

xDecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/xDecorr/", bc, "/", sys))))

GTwoDecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/GTwoDecorr/", bc, "/", sys))))

GFourDecorrData <- as.numeric(unlist(h5read("data.h5", paste0("/GFourDecorr/", bc, "/", sys))))

##### Thermalisation #####

# Some quick sanity checks

min(thermSweeps)  # Not too small (<1000)
max(thermSweeps)  # Not too big (>200000)

mean(accRateThermData) * 100 # Should be ~ 60 - 80%

# We need to do some work to the E0ThermData to get it into a form where we can plot the thermalisation curve.
# We know that the number of measurements during thermalisation per repeat is given by the number of sweeps divided by the thermalisation interval.

thermMeasures <- ceiling(thermSweeps / thermalisationInterval)  # Number of measurements during thermalisation per repeat
sum(thermMeasures) == length(E0ThermData) # Should return true

# We can split the E0ThermData into a list of vectors, where each vector corresponds to the thermalisation measurements for one repeat. 
E0Groups <- rep(seq_along(thermMeasures), thermMeasures)
E0Split <- split(E0ThermData, E0Groups)

# Then by finding the minimum length of these vectors, we can trim them all to the same length.
# This means we only plot the thermalisation curve up to the point where all repeats have measurements.
minLen <- min(sapply(E0Split, length))
E0Trimmed <- lapply(E0Split, `[`, 1:minLen)

# We then combine these trimmed vectors into a matrix, where each column corresponds to a repeat and each row corresponds to a measurement index.
E0Mat <- do.call(cbind, E0Trimmed)

# Then we take the average across repeats for each measurement index to get the average thermalisation curve.
E0Avg <- rowMeans(E0Mat)

ggplot(data.frame(index = 1:length(E0Avg), E0 = E0Avg), aes(x = index * thermalisationInterval, y = E0)) +
    geom_point() +
    labs(
      x = "Measurement index",
      y = "E0",
      title = paste("Thermalisation of", bc, sys)
    )

# Thermalisation of average position

#### FILL

##### Decorrelated data - Data we can use to find expected behaviour of the system and compare to analytical results #####

# Ground state energy

# Before we can work with decorrelated data, we need to split it into repeats. 
# We know that the number of measurements per repeat is given by measures.
E0Split <- split(E0DecorrData, rep(1:repeats, each = measures)) 

E0RepeatAvg <- sapply(E0Split, mean)  # Take mean per repeat

bins <- 7 # Change this to change the number of bins in the histogram

hE0 <- hist(E0RepeatAvg, breaks = bins, plot = FALSE) # Compute histogram breaks and bin width
binWidth <- diff(hE0$breaks)[1]

E0Norm <- hE0$counts / (sum(hE0$counts) * binWidth) # Normalised histogram counts (probability density)

continuousE0 <- seq(min(E0RepeatAvg), max(E0RepeatAvg), length.out = 200) # Continuous x values for normal curve
normDist <- dnorm(continuousE0, mean = mean(E0RepeatAvg), sd = sd(E0RepeatAvg))
dx <- diff(continuousE0)[1]
normDist <- normDist / sum(normDist * dx)  # Normalisation

# Histogram and normal curve
histPlot <- ggplot() +
  geom_histogram(aes(x = E0RepeatAvg, y = after_stat(density)), # after_stat(density) normalises the histogram to a density
                 bins = bins, fill = "skyblue", color = "black") +
  geom_line(aes(x = continuousE0, y = normDist),
            color = "red", linewidth = 1) +
  labs(title = "Ground State Energy Histogram",
       x = "Energy", y = "Probability Density") 

# QQ plot - this is a plot of the quantiles of our data against the quantiles of a normal distribution. 
# If the points lie approximately on a straight line, then our data is approximately normally distributed.
# Furthermore, we can perform a Shapiro-Wilk test for normality, which gives us a p-value. 
# If the p-value is above a certain threshold (usually 0.05), we fail to reject the null hypothesis that our
# data is normally distributed. (i.e. we can treat it as normally distributed for the purposes of our analysis).
shapiroTest <- shapiro.test(E0RepeatAvg)
qqPlot <- ggplot(data.frame(E0RepeatAvg), aes(sample = E0RepeatAvg)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = paste0("QQ Plot (Shapiro-Wilk p = ", round(shapiroTest$p.value, 4), ")"),
       x = "Theoretical Quantiles", y = "Sample Quantiles") 

# Combined plots side by side
grid.arrange(histPlot, qqPlot, ncol = 2)

histPlot # Show histogram and normal curve

qqPlot # Show QQ plot

mean(E0RepeatAvg) # Should be close to the expected ground state energy (0.5 for QHO, ~0.68 for DWP)

E0Range <- max(E0RepeatAvg) - min(E0RepeatAvg) 

monteCarloStandardError <- mean(E0RepeatAvg) * acceptableError

monteCarloStandardError * 2 # This is the 95% confidence interval for the mean, which should be smaller than the range of our data if we have a good estimate of the ground state energy.
E0Range > monteCarloStandardError * 2 # This should return true if we have a good estimate of the ground state energy.

mean(E0RepeatAvg) + monteCarloStandardError # Maximum E0 consistent with our data and error estimate
mean(E0RepeatAvg) - monteCarloStandardError # Minimum E0 consistent with our

# Wave function

bins <- 100

hist(xDecorrData, breaks = bins, main = "Histogram of positions",
     xlab = "position", ylab = "count")
h <- hist(xDecorrData, breaks = 100, plot = FALSE)

counts <- sum(h$counts)                             # Total counts which equals measures * lattice size
positionRange <- range(h$mids)[2] - range(h$mids)[1]  # Range of positions 
binWidth <- positionRange / bins                      # Width of each bin 
normFactor <- counts * binWidth                     # Normalization factor 
psi <- sqrt(h$counts / normFactor)                    # Normalized wave function 

if (sys == "QHO") {
  psiAnalytical <- exp(-(h$mids ^ 2) / 2)               # The analytical wavefunction of the QHO
} else if (sys == "AHO") {
  psiAnalytical <- exp(-(h$mids ^ 2) / 2)  # Approximation for the AHO to go here
} else if (sys == "DWP") {
  psiAnalytical <- exp(-(h$mids ^ 2) / 2)  # Put the WKB approximation of the DWP wavefunction here, atm it is just the QHO one
}

psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) # Normalising the analytical wavefunction 

ggplot(data.frame(x = h$mids, psi = psi), aes(x = x, y = psi)) +
  geom_line() +
  geom_point() +
  geom_line(aes(x = h$mids, y = psiAnalytical)) +
  labs(title = "Wave Function", x = "Position", y = "Psi")

ggplot(data.frame(x = h$mids, psi = psi), aes(x = x, y = psi)) +
  geom_line() +
  geom_point() +
  labs(title = "Wave Function", x = "Position", y = "Psi")


# Two point correlation function

expectedLength <- repeats * measures * pathLength
length(GTwoDecorrData) == expectedLength # Ensure xVec is the correct length

groupSize <- measures * pathLength

repeatIDs <- rep(1:repeats, each = groupSize)

GTwoDecorrDataSplit <- split(GTwoDecorrData, repeatIDs)

GTwoDecorrDataSplitMat <- lapply(GTwoDecorrDataSplit, function(v) {
  matrix(v, nrow = pathLength, ncol = measures)
})

# Compute the average correlation per repeat first
correlationList <- lapply(GTwoDecorrDataSplitMat, function(mat) {
  rowMeans(mat)  # Average over measurements for each lag
})

# Average over all repeats
correlationAvg <- Reduce("+", correlationList) / length(correlationList)

# Plot
dfCorr <- data.frame(
  lag = 0:(length(correlationAvg) - 1),
  correlation = correlationAvg
)

ggplot(dfCorr, aes(x = lag, y = correlation)) +
  geom_line(color = "#000000") +
  labs(
    title = paste("Two point decorrelation function for", bc, sys),
    x = "Time (index of the path)",
    y = "G_2(t, 0)"
  )


# E1 from the two-point correlator

successfulCounts <- 0
E1 <- 0

for (i in 1:10) {

  if (correlationAvg[i] <= 0 || correlationAvg[i + 1] <= 0) {
    message("Correlation function has non-positive values, cannot compute E1.")
  } 
  else {
    E1 <- E1 +
      mean(E0RepeatAvg) +
      log(correlationAvg[i] / correlationAvg[i + 1]) / latticeSpacing

    successfulCounts <- successfulCounts + 1
  }
}

E1 <- E1 / successfulCounts

E1

E1 - mean(E0Avg)

# Four point correlation function

expectedLength <- repeats * measures * pathLength
length(GFourDecorrData) == expectedLength # Ensure xVec is the correct length

groupSize <- measures * pathLength

repeatIDs <- rep(1:repeats, each = groupSize)

GFourDecorrDataSplit <- split(GFourDecorrData, repeatIDs)

GFourDecorrDataSplitMat <- lapply(GFourDecorrDataSplit, function(v) {
  matrix(v, nrow = pathLength, ncol = measures)
})

# Compute the average correlation per repeat first
correlationList <- lapply(GFourDecorrDataSplitMat, function(mat) {
  rowMeans(mat)  # Average over measurements for each lag
})

# Average over all repeats
correlationAvg <- Reduce("+", correlationList) / length(correlationList)

# Plot
dfCorr <- data.frame(
  lag = 0:(length(correlationAvg) - 1),
  correlation = correlationAvg
)

ggplot(dfCorr, aes(x = lag, y = correlation)) +
  geom_line(color = "#000000") +
  labs(
    title = paste("Four point decorrelation function for", bc, sys),
    x = "Time (index of the path)",
    y = "G_4(t, 0)"
  )


# E2 from the four point correlator

successfulCounts <- 0
E2 <- 0

for (i in 1:10) {

  if (correlationAvg[i] <= 0 || correlationAvg[i + 1] <= 0) {
    message("Correlation function has non-positive values, cannot compute E2.")
  } 
  else {
    E2 <- E2 +
      mean(E0RepeatAvg) +
      log(correlationAvg[i] / correlationAvg[i + 1]) / latticeSpacing

    successfulCounts <- successfulCounts + 1
  }
}

E2 <- E2 / successfulCounts

E2

E2 - mean(E0Avg)

# nolint end

