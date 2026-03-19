# nolint start

{
  library(ggplot2)
  library(gridExtra)  # side by side plots
  library(dplyr)
  library(rhdf5)
}

# Variables from the simulation, needed for interpreting the data
{
  # QHO variables
  m <- 1
  omega <- 1

  # AHO variables
  quarticFactor <- 1

  # DWP variables
  a <- 1.6
  lambda <- 1

  omegaDWP <- sqrt(8 * lambda * a^2)

  S_inst <- (2/3) * omegaDWP * (a^2) 

  alpha <- 1 / 12 # A complicated calculation performed in Zinn-Justin 1993 (or ABCs of Instantons) gives this value

  K <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5) # A prefactor for the splitting energy

  splittingEnergy <- K * exp(-S_inst)

  E0_inst <- 0.5 * omegaDWP - (splittingEnergy / 2)
  E1_inst <- 0.5 * omegaDWP + (splittingEnergy / 2)

  # Simulation values
  measures <- 50

  repeats <- 12

  pathLength <- 1600
  latticeSpacing <- 0.05
  beta <- pathLength * latticeSpacing

  thermalisationInterval <- 1
  acceptableError <- 0.01 
}

# Boundary conditions and system type
{
  bc <- "Periodic"
  # bc <- "Dirichlet"

  # sys <- "QHO"
  # sys <- "AHO"
  sys <- "DWP"
}

# Read data
{
  dataFile <- "data.h5"
    
  thermSweeps <- as.numeric(unlist(h5read(dataFile, paste0("/thermSweeps/", bc, "/", sys))))
  E0ThermData <- as.numeric(unlist(h5read(dataFile, paste0("/E0Therm/", bc, "/", sys))))
  accRateThermData <- as.numeric(unlist(h5read(dataFile, paste0("/accRateTherm/", bc, "/", sys))))
  E0Data <- as.numeric(unlist(h5read(dataFile, paste0("/E0/", bc, "/", sys))))
  accRateData <- as.numeric(unlist(h5read(dataFile, paste0("/accRate/", bc, "/", sys))))
  histogramData <- as.numeric(unlist(h5read(dataFile, paste0("/histogram/", bc, "/", sys))))
  instantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/instantons/", bc, "/", sys))))
  antiInstantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/antiInstantons/", bc, "/", sys))))
  GTwoData <- as.numeric(unlist(h5read(dataFile, paste0("/GTwo/", bc, "/", sys))))
  GFourData <- as.numeric(unlist(h5read(dataFile, paste0("/GFour/", bc, "/", sys))))
}

##### Thermalisation #####

min(thermSweeps)
max(thermSweeps)

# Acceptance rate
mean(accRateThermData) * 100 # Should be ~ 50 - 80%

# Getting average thermalisation of the ground state energy

{
  # Compute differences between consecutive points
  diffE <- diff(E0ThermData)

  # Detect new repeat starts by large negative jumps in energy
  threshold <- -0.3  # adjust for different systems
  new_repeat_indices <- which(diffE < threshold) + 1  # +1 because diff shifts index

  # Include the very first point as start of first repeat
  repeat_starts <- c(1, new_repeat_indices)

  # Assign repeat IDs to each data point
  repeat_id <- rep(seq_along(repeat_starts), 
                  times = c(diff(c(repeat_starts, length(E0ThermData)+1))))

  # Build a data frame for plotting
  df_detected <- data.frame(
    sweep_index = seq_along(E0ThermData),
    E0 = E0ThermData,
    rep = factor(repeat_id)
  )

  # Compute per-repeat lengths
  repeat_lengths <- diff(c(repeat_starts, length(E0ThermData) + 1))

  # Find the minimum length to align all repeats
  min_len <- min(repeat_lengths)

  # Extract and trim each repeat
  E0Trimmed <- lapply(seq_along(repeat_starts), function(i) {
    start <- repeat_starts[i]
    end <- start + min_len - 1
    E0ThermData[start:end]
  })

  # Combine into a matrix: rows = sweeps, columns = repeats
  E0Mat <- do.call(cbind, E0Trimmed)

  # Compute average thermalisation curve across repeats
  E0Avg <- rowMeans(E0Mat)

  # Build a sweep index for plotting 
  sweep_index <- seq_len(min_len)  # simple sweep numbers
}

# Plot average thermalisation
ggplot(data.frame(sweep = sweep_index, E0 = E0Avg), aes(x = sweep, y = E0)) +
  geom_line(color = "blue", size = 1) +
  labs(
    x = "MC Sweep (aligned)",
    y = "Average E0",
    title = "Average Thermalisation Across Repeats"
  ) +
  theme_minimal()

##### Decorrelated data - Data we can compare against theory #####

### Ground state energy

# Histogram and shapiro test
{
  bins <- 7 

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

qqPlot # Show QQ plot

# Calculating E0 and its error
{
  E0 <- mean(E0RepeatAvg) # Should be close to the expected ground state energy (0.5 for QHO, ~0.68 for DWP)

  E0Range <- max(E0RepeatAvg) - min(E0RepeatAvg)
  monteCarloStandardError <- mean(E0RepeatAvg) * acceptableError

  monteCarloStandardError * 2 # This is the 95% confidence interval for the mean, which should be smaller than the range of our data if we have a good estimate of the ground state energy.
  E0Range > monteCarloStandardError * 2 # This should return true if we have a good estimate of the ground state energy.
}
E0; mean(E0RepeatAvg) + monteCarloStandardError; mean(E0RepeatAvg) - monteCarloStandardError 

### Wave function

# Histogram data frame creation
{
  numBins <- 100 # Same as in C++

  # Histogram range
  if (sys == "FP") {
    xMax <- 5.0
    xMin <- -5.0
  } else if (sys == "QHO" || sys == "AHO") {
    sigmaQHO <- 1 / sqrt(2 * m * omega)
    xMax <- ceiling(4 * sigmaQHO) + 1  # 4 standard deviations of analytic ground state plus padding
    xMin <- -xMax
  } else if (sys == "DWP") {
    sigmaDWP <- 1 / sqrt(omegaDWP)
    xMax <- ceiling(a + 4 * sigmaDWP) + 1
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
    psiAnalytical <- exp(-(omegaDWP / 2) * (x_values + a)^2) + exp(-(omegaDWP / 2) * (x_values - a)^2) # Approximate by 2 Gaussians
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

  LHS <- 1; RHS <- 800

  correlatorRatios <- numeric(RHS - LHS + 1) # To store the log ratios for each lag

  for (i in LHS:RHS) {
    if (GTwoData[i] <= 0 || GTwoData[i + 1] <= 0) {
      message("Correlation function has non-positive values, cannot compute E1.")
    } 
    else {
      E1 <- E1 + mean(E0RepeatAvg) + log(GTwoData[i] / GTwoData[i + 1]) / latticeSpacing
      correlatorRatios[i - LHS + 1] <- log(GTwoData[i] / GTwoData[i + 1]) / latticeSpacing
      successfulCounts <- successfulCounts + 1
    }
  }
  E1 <- E1 / successfulCounts; E1
}

ggplot(data.frame(lag = 1:(length(correlatorRatios)), ratio = correlatorRatios), aes(x = lag, y = ratio)) +
  geom_line(color = "blue") +
  labs(title = "Log Ratio of Two Point Correlation Function",
       x = "Lag", y = "log(G_2(t) / G_2(t+1))")

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
beta / mean(instantonsData)

1 / omegaDWP

# Condition for good instanton sampling

exp(-S_inst) # Must be << 1

# Some extra calculations 

E0_inst
E1_inst
E0
E1
splittingEnergy
E1 - E0

# nolint end

3.63^2
