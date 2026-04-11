# nolint start

{
  library(ggplot2)
  library(gridExtra)  # side by side plots
  library(dplyr)
  library(rhdf5)
}

# Boundary conditions and system type
{
  bc <- "Periodic"
  # bc <- "Dirichlet"

  sys <- "QHO"
  # sys <- "AHO"
  # sys <- "DWP"
}

# Read data
{
  dataFile <- "data.h5"
  
  E0ThermData <- as.numeric(unlist(h5read(dataFile, paste0("/E0Therm/", bc, "/", sys))))
  accRateThermData <- as.numeric(unlist(h5read(dataFile, paste0("/accRateTherm/", bc, "/", sys))))
  E0Data <- as.numeric(unlist(h5read(dataFile, paste0("/E0/", bc, "/", sys))))
  accRateData <- as.numeric(unlist(h5read(dataFile, paste0("/accRate/", bc, "/", sys))))
  histogramData <- as.numeric(unlist(h5read(dataFile, paste0("/histogram/", bc, "/", sys))))
  instantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/instantons/", bc, "/", sys))))
  antiInstantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/antiInstantons/", bc, "/", sys))))
  GTwoData <- as.numeric(unlist(h5read(dataFile, paste0("/GTwo/", bc, "/", sys))))
  GFourData <- as.numeric(unlist(h5read(dataFile, paste0("/GFour/", bc, "/", sys))))
  headerData <- as.numeric(unlist(h5read(dataFile, paste0("/headerInfo/", bc, "/", sys))))
}

# Variables from the simulation
{
  pathLength <- headerData[1]
  latticeSpacing <- headerData[2]
  epsilon <- headerData[3]
  accRateInterval <- headerData[4]
  decorrSweeps <- headerData[5]
  thermSweeps <- headerData[6]
  thermInterval <- headerData[7]
  measures <- headerData[8]
  repeats <- headerData[9]
  numBins <- headerData[10]
  
  beta <- pathLength * latticeSpacing
  thermMeasures <- thermSweeps / thermInterval

  if (sys == "QHO") {
    mQHO <- headerData[11]
    omegaQHO <- headerData [12]
  } else if (sys == "AHO") {
    lambdaAHO <- headerData[11] # Should use more parameters for the AHO
  } else if (sys == "DWP") {
    wellCentres <- headerData[11]
    lambdaDWP <- headerData [12]
    omegaDWP <- sqrt(8 * lambdaDWP * wellCentres^2)
    S_inst <- (2/3) * omegaDWP * (wellCentres^2) 
    alpha <- 1 / 12 # A complicated calculation performed in Zinn-Justin 1993 (or ABCs of Instantons) gives this value
    K <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5) # A prefactor for the splitting energy
    splittingEnergy <- K * exp(-S_inst)
    E0_inst <- 0.5 * omegaDWP - (splittingEnergy / 2)
    E1_inst <- 0.5 * omegaDWP + (splittingEnergy / 2)
  }
}

### Thermalisation ###

# Acceptance rate
mean(accRateThermData) * 100 # Should be between 50 and 80%

# Getting average thermalisation of the ground state energy
{
  # Split into repeats, each column is one repeat
  E0Mat <- matrix(E0ThermData, nrow = thermMeasures, ncol = repeats, byrow = FALSE)

  # Compute average thermalisation across repeats
  E0ThermAvgs <- rowMeans(E0Mat)

  # Sweep index for plotting
  sweepIndex <- seq_len(thermMeasures)

  # Data frame for plotting
  E0ThermDF <- data.frame(sweep = sweepIndex, E0 = E0ThermAvgs)
}

# Plot thermalisation
ggplot(data.frame(sweep = sweepIndex, E0 = E0ThermAvgs), aes(x = sweep * thermInterval, y = E0)) +
  geom_line(color = "blue", size = 1) +
  labs(
    x = "Measure index",
    y = "Average Ground state energy",
    title = "Average Thermalisation Across Measures"
  ) +
  theme_minimal()
# Expect the ground state energy to taper off for a positive test of whether thermalisation has been accomplished

### Ground state energy ###

# Decorrelation
{
  # if (measures > 50) {
  #   decorrCheck <- 50
  # } else {
  #   decorrCheck <- measures
  # }
  decorrCheck <- measures
  measureIndex <- seq_len(decorrCheck)
  ggplot(data.frame(sweep = measureIndex, E0 = E0Data[1:decorrCheck]), aes(x = sweep, y = E0)) + 
    geom_line(color = "blue", size = 1) +
    labs(
      x = "MC Sweep",
      y = "Average E0",
      title = "Average Thermalisation Across Repeats"
    ) +
    theme_minimal()
  # For decorrelated data, look for low correlation between measure index and average ground state energy
}

# Further proof of decorrelation
{
  acfResult <- acf(E0Data[1:measures], lag.max = measures - 1, plot = TRUE)
  acfVals <- acfResult$acf
  mean(abs(as.vector(acfVals)[2:length(acfVals)])) # For decorrelated data, we want low acf values (roughly < 0.1)
}

# Histogram and shapiro test
{
  bins <- 7   # Can be increased for high number of repeats

  E0Split <- split(E0Data, rep(1:repeats, each = measures)) # Split E0 into repeats
  E0RepeatAvg <- sapply(E0Split, mean)  # Take mean per repeat

  hE0 <- hist(E0RepeatAvg, breaks = bins, plot = FALSE) # Compute histogram breaks and bin width
  binWidth <- diff(hE0$breaks)[1]

  E0Norm <- hE0$counts / (sum(hE0$counts) * binWidth) # Normalised histogram counts (probability density)

  continuousE0 <- seq(min(E0RepeatAvg), max(E0RepeatAvg), length.out = 200) # Continuous x values for normal curve
  normDist <- dnorm(continuousE0, mean = mean(E0RepeatAvg), sd = sd(E0RepeatAvg))
  dx <- diff(continuousE0)[1]
  normDist <- normDist / sum(normDist * dx)  # Normalisation

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
}

grid.arrange(histPlot, qqPlot, ncol = 2)  # Combined plots side by side

histPlot # Show histogram and normal curve

qqPlot # Show QQ plot (We want p > 0.05)

# Calculating E0 and its error
{
  E0 <- mean(E0RepeatAvg) # Should be close to the expected ground state energy (0.5 for QHO, ~0.68 for DWP)

  E0StandardError <- sd(E0RepeatAvg) / sqrt(length(E0RepeatAvg))  # Standard error
}
E0; mean(E0RepeatAvg) + E0StandardError; mean(E0RepeatAvg) - E0StandardError 
((E0 - 0.5) / 0.5) * 100
### Wave function ###

# Histogram data frame creation
{
  # Histogram range
  if (sys == "QHO") {
    sigmaQHO <- 1 / sqrt(2 * mQHO * omegaQHO)
    xMax <- ceiling(4 * sigmaQHO) + 1  # 4 standard deviations of analytic ground state plus padding
    xMin <- -xMax
  } else if (sys == "AHO") {
    sigmaAHO <- 1 / sqrt(omegaAHO)
    xMax <- ceiling(4 * sigmaAHO) + 1
    xMin <- -xMax
  } else if (sys == "DWP") {
    sigmaDWP <- 1 / sqrt(omegaDWP)
    xMax <- ceiling(wellCentres + 4 * sigmaDWP) + 1
    xMin <- -xMax
  }

  binWidth <- (xMax - xMin) / numBins
  x_values <- seq(xMin + binWidth/2, xMax - binWidth/2, length.out = numBins)

  hist_matrix <- matrix(histogramData, nrow = repeats, ncol = numBins, byrow = TRUE)
  hist_avg <- colMeans(hist_matrix)

  prob_density <- hist_avg / (sum(hist_avg) * binWidth)

  hist_df <- data.frame(
    x = x_values,
    probability = prob_density
  )
}

# Plot histogram
ggplot(hist_df, aes(x = x, y = probability)) +
  geom_col(width = binWidth, fill = "steelblue") +   # bar plot
  # Or use geom_line() + geom_point() for a line plot:
  # geom_line(color = "steelblue", size = 1) +
  # geom_point(color = "darkblue", size = 2) +
  labs(title = paste(sys, "position histogram"),
       x = "x",
       y = "Probability density") +
  theme_minimal(base_size = 14)

psi <- sqrt(prob_density)

# Finding the wavefunction from the histogram
{
  wave_df <- data.frame(
    x = x_values,
    psi = psi
  )

  # Overlay analytical wavefunction 
  if (sys == "QHO") {
    psiAnalytical <- exp(-(x_values^2) / 2)
  } else if (sys == "AHO") {
    psiAnalytical <- exp(-(x_values^2) / 2)  # Placeholder
  } else if (sys == "DWP") {
    psiAnalytical <- exp(-(omegaDWP / 2) * (x_values + wellCentres)^2) + exp(-(omegaDWP / 2) * (x_values - wellCentres)^2) # Approximate by 2 Gaussians
  }

  # Normalise the analytical wavefunction 
  psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) 

  wave_df$psiAnalytical <- psiAnalytical
}

# Plot the wavefunction
ggplot(wave_df, aes(x = x)) +
  geom_line(aes(y = psi), color = "blue") +
  geom_point(aes(y = psi), color = "darkblue", size = 1.5) +
  geom_line(aes(y = psiAnalytical), color = "red", linetype = "dashed") +
  labs(title = paste(sys, "Wave Function"),
       x = "Position", y = "Psi") +
  theme_minimal(base_size = 14)

### Two point correlation function

dfCorr <- data.frame(
  lag = 0:(length(GTwoData) - 1),
  correlation = GTwoData
)
sqrt(min(GTwoData))
# GTwoData

ggplot(dfCorr, aes(x = lag, y = correlation)) +
  geom_line(color = "#000000") +
  labs(
    title = paste("Two point decorrelation function for", bc, sys),
    x = "Time (index of the path)",
    y = "G_2(t, 0)"
  )

# Four point correlation function

dfCorr <- data.frame(
  lag = 0:(length(GFourData) - 1),
  correlation = GFourData
)

ggplot(dfCorr, aes(x = lag, y = correlation)) +
  geom_line(color = "#000000") +
  labs(
    title = paste("Four point decorrelation function for", bc, sys),
    x = "Time (index of the path)",
    y = "G_4(t, 0)"
  )

### Excited energy states

# E1

{
  successfulCounts <- 0; E1 <- 0

  LHS <- 1; RHS <- pathLength / 2

  correlatorRatios <- numeric(RHS - LHS + 1) # To store the log ratios for each lag

  for (i in LHS:RHS) {
    if (GTwoData[i] <= 0 || GTwoData[i + 1] <= 0) {
      message("Correlation function has non-positive values, cannot compute logarithmic ratio.")
    } 
    else {
      E1 <- E1 + mean(E0RepeatAvg) + log(GTwoData[i] / GTwoData[i + 1]) / latticeSpacing
      correlatorRatios[i - LHS + 1] <- log(GTwoData[i] / GTwoData[i + 1]) / latticeSpacing
      successfulCounts <- successfulCounts + 1
    }
  }
  E1 <- E1 / successfulCounts; E1
}

ggplot(data.frame(lag = 1:(length(correlatorRatios)), ratio = correlatorRatios), aes(x = lag * latticeSpacing, y = ratio)) +
  geom_line(color = "blue") +
  labs(title = "Log Ratio of Two Point Correlation Function",
       x = "Lag", y = "log(G_2(t) / G_2(t+1))")

find_flat_region <- function(x, y, window = 10) {
  n <- length(x)
  slopes <- rep(NA, n - window + 1)
  
  for (i in 1:(n - window + 1)) {
    fit <- lm(y[i:(i+window-1)] ~ x[i:(i+window-1)])
    slopes[i] <- abs(coef(fit)[2])
  }
  
  best_idx <- which.min(slopes)
  
  list(
    start = best_idx,
    end = best_idx + window - 1,
    slope = slopes[best_idx]
  )
}

find_flat_region(seq(1:length(correlatorRatios)), correlatorRatios)

min(correlatorRatios)
which.min(correlatorRatios)

E1 - E0

{
  successfulCounts <- 0; E2 <- 0

  for (i in 1:10) {
    if (GFourData[i] <= 0 || GFourData[i + 1] <= 0) {
      message("Correlation function has non-positive values, cannot compute E2.")
    } 
    else {
      E2 <- E2 + mean(E0RepeatAvg) + log(GFourData[i] / GFourData[i + 1]) / latticeSpacing
      successfulCounts <- successfulCounts + 1
    }
  }
  E2 <- E2 / successfulCounts; E2
}

E2 - E0

# Tunnelling data

head(instantonsData)
head(antiInstantonsData)
mean(instantonsData)
mean(instantonsData) / beta

exp(- S_inst)

1 / omegaDWP

# Condition for good instanton sampling

exp(-S_inst) # Must be << 1

S_inst_approximation <- log(pathLength / mean(instantonsData)) # This is an approximation for the instanton action based on the number of tunnelling events observed in the simulation (N is the total number of paths sampled, tunnelingEvents is the number of paths that exhibited tunnelling)
S_inst_approximation * omegaDWP / pi * sqrt(S_inst_approximation)
S_inst

# Some extra calculations 

E0_inst
E1_inst
E0
E1
splittingEnergy
E1 - E0

# nolint end