# nolint start

{
  library(ggplot2)
  library(gridExtra)  # side by side plots
  library(dplyr)
  library(rhdf5)
  library(minpack.lm) # better exponential fitting
}

# Boundary conditions and system type
{
  bc <- "Periodic"
  # bc <- "Dirichlet"

  # sys <- "QHO"
  sys <- "DWP"
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

  diagEnergiesData <- read.csv("DWP diagonalisation/DWP diagonalisation/energies.csv")
  diagWFData <- read.csv("DWP diagonalisation/DWP diagonalisation/wavefunctions.csv")
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
    omegaQHO <- headerData[12]
  } else if (sys == "DWP") {
    wellCentres <- headerData[11]
    lambdaDWP <- headerData[12]

    # Diagonalisation Energies
    energy_row <- as.integer(round((wellCentres - 0.9) * 10))
    E0_diag <- diagEnergiesData[energy_row, 2]
    E1_diag <- diagEnergiesData[energy_row, 3]
    Split_diag <- diagEnergiesData[energy_row, 4]

    # WKB analysis
    omegaDWP <- sqrt(8 * lambdaDWP * wellCentres^2)
    S_inst <- (2/3) * omegaDWP * (wellCentres^2) 
    alpha <- 1 / 12 # A complicated calculation performed in Zinn-Justin 1993 (or ABCs of Instantons) gives this value
    K <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5) # A prefactor for the splitting energy
    Split_WKB <- K * exp(-S_inst)
    E0_WKB <- 0.5 * omegaDWP - (Split_WKB / 2)
    E1_WKB <- 0.5 * omegaDWP + (Split_WKB / 2)
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
{
  if (sys == "QHO") {
    ggplot(data.frame(sweep = sweepIndex, E0 = E0ThermAvgs), aes(x = sweep * thermInterval, y = E0)) +
      geom_line(color = "blue", linewidth = 0.6) +
      labs(x = "Sweep", y = expression("Ground State Energy" ~ E[0]),
        title = "Average Thermalisation of the QHO") +
      theme_minimal(base_size = 14)
    ggsave("Figures/QHO_thermalisation.png", width = 6, height = 4, dpi = 1200)
  } else if (sys == "DWP") {
    ggplot(data.frame(sweep = sweepIndex, E0 = E0ThermAvgs), aes(x = sweep * thermInterval, y = E0)) +
      geom_line(color = "blue", linewidth = 0.6) +
      labs(x = "Sweep", y = expression("Ground State Energy" ~ E[0]),
        title = "Average Thermalisation of the DWP") +
      theme_minimal(base_size = 14)
    ggsave("Figures/DWP_thermalisation.png", width = 6, height = 4, dpi = 1200)
  }
}

### Ground state energy ###

# Decorrelation
{
  if (measures > 50) {
    decorrCheck <- 50
  } else {
    decorrCheck <- measures
  }
  measureIndex <- seq_len(decorrCheck)
  if (sys == "QHO") {
    ggplot(data.frame(sweep = measureIndex, E0 = E0Data[1:decorrCheck]), aes(x = sweep, y = E0)) + 
      geom_line(color = "blue", size = 0.6) +
      labs(x = "Measure", y = expression("Ground State Energy" ~ E[0]), 
        title = "Average Decorrelation of the QHO") +
      theme_minimal(base_size = 14)
    ggsave("Figures/QHO_decorrelation.png", width = 6, height = 4, dpi = 1200)  
  } else if (sys == "DWP") {
    ggplot(data.frame(sweep = measureIndex, E0 = E0Data[1:decorrCheck]), aes(x = sweep, y = E0)) + 
      geom_line(color = "blue", size = 0.6) +
      labs(x = "Measure", y = expression("Ground State Energy" ~ E[0]), 
        title = "Average Decorrelation of the QHO") +
      theme_minimal(base_size = 14)
    ggsave("Figures/DWP_decorrelation.png", width = 6, height = 4, dpi = 1200)  
  }
}

# Further proof of decorrelation - auto correlation function
{
  acfResult <- acf(E0Data[1:measures], lag.max = measures - 1, plot = TRUE)
  acfVals <- acfResult$acf
  mean(abs(as.vector(acfVals)[2:length(acfVals)])) # For decorrelated data, we want low acf values (roughly < 0.1)
}

# Save acf plot
{
  if (sys == "QHO")  {
    png("Figures/acf_plot_QHO.png", width = 800, height = 600)

    acfResult <- acf(E0Data[1:measures], lag.max = measures - 1, plot = FALSE)

    plot(acfResult,
        main = "Autocorrelation of Ground State Energy for the QHO",
        xlab = "Lag (sweeps)",
        ylab = "Autocorrelation")

    dev.off()
    } else if (sys == "DWP") {
      png("Figures/acf_plot_DWP.png", width = 800, height = 600)

      acfResult <- acf(E0Data[1:measures], lag.max = measures - 1, plot = FALSE)

      plot(acfResult,
          main = "Autocorrelation of Ground State Energy for the DWP",
          xlab = "Lag (sweeps)",
          ylab = "Autocorrelation")

      dev.off()
    }
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

  if (sys == "QHO") {
    histPlot <- ggplot() +
      geom_histogram(aes(x = E0RepeatAvg, y = after_stat(density)), # after_stat(density) normalises the histogram to a density
                    bins = bins, fill = "skyblue", color = "black") +
      geom_line(aes(x = continuousE0, y = normDist),
                color = "red", linewidth = 1) +
      labs(title = "Ground State Energy Histogram for the QHO",
          x = expression("Ground State Energy " ~ E[0]), y = "Probability Density") 

    shapiroTest <- shapiro.test(E0RepeatAvg)
    qqPlot <- ggplot(data.frame(E0RepeatAvg), aes(sample = E0RepeatAvg)) +
      stat_qq() +
      stat_qq_line(color = "red") +
      labs(title = paste0("QQ Plot (Shapiro-Wilk p = ", round(shapiroTest$p.value, 4), ")"),
          x = "Theoretical Quantiles", y = "Sample Quantiles")
    png("Figures/QHO_Histogram_qq.png", width = 1200, height = 800)
    grid.arrange(histPlot, qqPlot, ncol = 2)
    dev.off()
  }
  else if (sys == "DWP") {
    histPlot <- ggplot() +
      geom_histogram(aes(x = E0RepeatAvg, y = after_stat(density)), # after_stat(density) normalises the histogram to a density
                    bins = bins, fill = "skyblue", color = "black") +
      geom_line(aes(x = continuousE0, y = normDist),
                color = "red", linewidth = 1) +
      labs(title = "Ground State Energy Histogram for the DWP",
          x = expression("Ground State Energy " ~ E[0]), y = "Probability Density") 

    shapiroTest <- shapiro.test(E0RepeatAvg)
    qqPlot <- ggplot(data.frame(E0RepeatAvg), aes(sample = E0RepeatAvg)) +
      stat_qq() +
      stat_qq_line(color = "red") +
      labs(title = paste0("QQ Plot (Shapiro-Wilk p = ", round(shapiroTest$p.value, 4), ")"),
          x = "Theoretical Quantiles", y = "Sample Quantiles")
    png("Figures/DWP_Histogram_qq.png", width = 1200, height = 800)
    grid.arrange(histPlot, qqPlot, ncol = 2)
    dev.off()
  }
}

# grid.arrange(histPlot, qqPlot, ncol = 2)  # Combined plots side by side

# histPlot # Show histogram and normal curve

# qqPlot # Show QQ plot

# Calculating E0 and its error
{
  E0 <- mean(E0RepeatAvg) # Should be close to the expected ground state energy (0.5 for QHO, ~0.68 for DWP)

  E0StandardError <- sd(E0RepeatAvg) / sqrt(length(E0RepeatAvg))  # Standard error
}

E0; mean(E0RepeatAvg) + E0StandardError; mean(E0RepeatAvg) - E0StandardError 

# Percent error in E0
{
  if (sys == "QHO") {
    E0_error <- abs(((E0 - 0.5) / 0.5) * 100) 
  }
  else if (sys == "DWP") {
    E0_error <- abs(((E0 - E0_diag) / E0_diag) * 100) # Approximate error in E0 by the most accurate method known
  }
  E0_error
}

### Wave function ###

# Histogram data frame creation
{
  # Histogram range
  if (sys == "QHO") {
    sigmaQHO <- 1 / sqrt(2 * mQHO * omegaQHO)
    xMax <- ceiling(4 * sigmaQHO) + 1  # 4 standard deviations of analytic ground state plus padding
    xMin <- -xMax
  }
  else if (sys == "DWP") {
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

# Finding the wave function from the histogram
{
  psi <- sqrt(prob_density)

  wave_df <- data.frame(
    x = x_values,
    psi = psi
  )

  # Overlay analytical wave function 
  if (sys == "QHO") {
    psiAnalytical <- exp(-(x_values^2) / 2)
  } else if (sys == "DWP") {
    psiAnalytical <- exp(-(omegaDWP / 2) * (x_values + wellCentres)^2) + exp(-(omegaDWP / 2) * (x_values - wellCentres)^2) # Approximate by 2 Gaussians
  }

  # Normalise the analytical wave function 
  psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) 

  wave_df$psiAnalytical <- psiAnalytical
  
  ## Diagonalisation wave function
  if (sys == "DWP")  { 
    f_target <- wellCentres
    diagWFData <- diagWFData %>%
      filter(abs(f - f_target) < 1e-6)
  }
}

# Comparing wave functions

{
  if (sys == "QHO") {
    ggplot(wave_df, aes(x = x)) +
      # MCMC wave function
      geom_point(aes(y = psi, color = "MCMC"), size = 1) +

      # Exact wave function
      geom_line(aes(y = psiAnalytical, color = "Exact")) +

      labs(title = paste(sys, "Quantum Harmonic Oscillator Ground State Wave Function"), x = "Position", y = "Psi") +

      scale_y_continuous(
        name = expression("Wavefunction " ~ psi[0](x))) +

      scale_color_manual(
        name = "",
        values = c("MCMC" = "red", "Exact" = "black")) + 

      theme_minimal(base_size = 14)
    ggsave("Figures/QHO_wave_functions.png", width = 8, height = 4, dpi = 1200)  
  } else if (sys == "DWP") {
    ggplot(wave_df, aes(x = x)) +
      # MCMC wave function
      geom_point(aes(y = psi, color = "MCMC"), size = 1) +

      # Diagonalisation wave function
      geom_line(data = diagWFData,
                aes(x = x, y = psi0, color = "Diagonalised"),
                linewidth = 0.7) +

      # WKB approximation
      geom_line(aes(y = psiAnalytical, color = "WKB"), linetype = "dashed") +
      
      labs(title = "Double-Well Potential Ground State Wave Function", x = "Position x") +

      # coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 0.8)) + 

      scale_y_continuous(
        name = expression("Wavefunction " ~ psi[0](x))) +

      scale_color_manual(
        name = "",
        values = c("MCMC" = "red", "Diagonalised" = "black", "WKB" = "blue")) + 
  
      theme_minimal(base_size = 14)
    ggsave("Figures/DWP_wave_functions.png", width = 8, height = 4, dpi = 1200)
  }
}

### Two point correlation function

# Exponential fit

{
  noiseless_region_2 <- 200 # Adjust to the length of the noiseless region
  dfCorr <- data.frame(
    lag = 0:(noiseless_region_2 - 1),
    correlation = GTwoData[1:noiseless_region_2]
  )
  fit_data <- subset(dfCorr, lag >= 0 & lag <= noiseless_region_2)
  if (sys == "QHO") {
    DeltaE_guess <- 1.0 # Hard coded energy splitting for QHO with m = omega = 1
  } else if (sys == "DWP") {
    DeltaE_guess <- Split_diag
  }
  A_guess <- 0.098
  c_guess <- min(fit_data$correlation)
  
  fit <- nlsLM(
    correlation ~ A * exp(-DeltaE * lag * latticeSpacing) + c,
    data = fit_data,
    start = list(A = A_guess, c= c_guess, DeltaE = DeltaE_guess),
    control = nls.lm.control(maxiter = 1024)
  )
  dfCorr$fit <- predict(fit, newdata = dfCorr)
  coef(fit)
}

{
  if (sys == "QHO") {
    ggplot(dfCorr, aes(x = lag * latticeSpacing)) +
      geom_point(aes(y = correlation, color = "MCMC Correlator"), size = 0.5) +
      geom_line(aes(y = fit, color = "Exponential fit")) +
      labs(title = "QHO Two-Point Correlator with Exponential fit", y = "Two-Point Correlator G(t)", x = "Lag t") +
      scale_color_manual(
        name = "",
        values = c("MCMC Correlator" = "black", "Exponential fit" = "red"))
    ggsave("Figures/QHO_two_point_correlator.png", width = 7, height = 4, dpi = 1200)
  } else if (sys == "DWP") {
    ggplot(dfCorr, aes(x = lag * latticeSpacing)) +
      geom_point(aes(y = correlation, color = "MCMC Correlator"), size = 0.5) +
      geom_line(aes(y = fit, color = "Exponential fit")) +
      labs(title = "DQP Two-Point Correlator with Exponential fit", y = "Two-Point Correlator G(t)", x = "Lag t") +
      scale_color_manual(
        name = "",
        values = c("MCMC Correlator" = "black", "Exponential fit" = "red"))
    ggsave("Figures/DWP_two_point_correlator.png", width = 7, height = 4, dpi = 1200)
  }
}

# Four point correlation function
{
  noiseless_region_4 <- 40 # Adjust to the length of the noiseless region
  dfCorr <- data.frame(
    lag = 0:(noiseless_region_4 - 1),
    correlation = GFourData[0:noiseless_region_4]
  )
}

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

  LHS <- 1; RHS <- noiseless_region_2

  correlatorRatios <- numeric(RHS - LHS + 1) # To store the log ratios for each lag

  for (i in LHS:RHS) {
    if (GTwoData[i] <= 0 || GTwoData[i + 1] <= 0) {
      message("Correlation function has non-positive values, cannot compute log ratio.")
    } 
    else {
      E1 <- E1 + E0 + log(GTwoData[i] / GTwoData[i + 1]) / latticeSpacing
      correlatorRatios[i - LHS + 1] <- log(GTwoData[i] / GTwoData[i + 1]) / latticeSpacing
      successfulCounts <- successfulCounts + 1
    }
  }
  E1 <- E1 / successfulCounts; E1
}

{
  if (sys == "QHO") {   
    ggplot(data.frame(lag = 1:(length(correlatorRatios)), ratio = correlatorRatios), aes(x = lag * latticeSpacing, y = ratio)) +
      geom_line(aes(color = "Log Ratio")) +
      geom_hline(data = data.frame(y = E1 - E0), aes(yintercept = y, color = "Splitting Energy"), linetype = "dashed") +
      labs(title = "Finding the splitting energy for the QHO",
          x = "Lag t", y = expression("Log(" ~ G[2](t) ~ ") / Log(" ~ G[2](t+a) ~ ")")) +
      scale_color_manual(name = "", values = c("Log Ratio" = "blue", "Splitting Energy" = "red"))
    ggsave("Figures/QHO_log_ratio.png", width = 6, height = 4, dpi = 1200)

  } else if (sys == "DWP") {
    ggplot(data.frame(lag = 1:(length(correlatorRatios)), ratio = correlatorRatios), aes(x = lag * latticeSpacing, y = ratio)) +
      geom_line(aes(color = "Log Ratio")) +
      geom_hline(data = data.frame(y = E1 - E0), aes(yintercept = y, color = "Splitting Energy"), linetype = "dashed") +
      labs(title = "Finding the splitting energy for the DWP",
          x = "Lag t", y = expression("Log(" ~ G[2](t) ~ ") / Log(" ~ G[2](t+a) ~ ")")) +
      scale_color_manual(name = "", values = c("Log Ratio" = "blue", "Splitting Energy" = "red"))
    ggsave("Figures/DWP_log_ratio.png", width = 6, height = 4, dpi = 1200)
  } 
}

E1 - E0
Split_diag

{
  successfulCounts <- 0; E2 <- 0

  for (i in 1:10) {
    if (GFourData[i] <= 0 || GFourData[i + 1] <= 0) {
      message("Correlation function has non-positive values, cannot compute log ratio.")
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

mean(instantonsData)

# Conditions for good instanton sampling

exp(-S_inst) # Must be << 1
beta * Split_diag # Should be > 10

# An approximation for the instanton action based on the number of tunnelling events  
{
  S_inst_approximation <- log(pathLength / mean(instantonsData)) 
  S_inst_approximation <- S_inst_approximation * omegaDWP / pi * sqrt(S_inst_approximation)
} 

S_inst; S_inst_approximation

# Comparing energies against diagonalisation method and WKB

E0; E0_diag; E0_WKB

E1; E1_diag; E1_WKB

E1 - E0; Split_diag; Split_WKB 

##### Discretisation error analysis (DWP)

{
  EA_data <- read.csv("Discretisation error DWP.csv")
  EA_data <- EA_data %>%
    filter(!is.na(E0MCMC), !is.na(E1MCMC))

  theory_df <- EA_data %>%
    distinct(Separation, E0Diag, E0WKB, E1Diag, E1WKB, SplitDiag, SplitWKB) # CamelCase for vals in csv
}

ggplot(EA_data, aes(x = 1/(a * a))) +

  # MCMC
  geom_point(aes(y = E0MCMC, color = factor(Separation))) +
  geom_line(aes(y = E0MCMC, color = factor(Separation))) +

  # Diagonalisation horizontal lines
  geom_hline(data = theory_df,
             aes(yintercept = E0Diag, linetype = "Diagonalisation"),
             color = "red") +

  # WKB horizontal lines
  # geom_hline(data = theory_df,
  #            aes(yintercept = E0WKB, linetype = "WKB"),
  #            color = "black") +

  scale_linetype_manual(values = c(
    "Diagonalisation" = "dashed",
    "WKB" = "dotted"
  )) +

  labs(
    title = "E0 vs 1/(a^2)",
    linetype = "Theory"
  )

ggplot(EA_data, aes(x = 1/(a * a))) +

  # MCMC
  geom_point(aes(y = SplitMCMC, color = factor(Separation))) +
  geom_line(aes(y = SplitMCMC, color = factor(Separation))) +

  # Diagonalisation horizontal lines
  geom_hline(data = theory_df,
             aes(yintercept = SplitDiag, linetype = "Diagonalisation"),
             color = "red") +

  # WKB horizontal lines
  geom_hline(data = theory_df,
             aes(yintercept = SplitWKB, linetype = "WKB"),
             color = "black") +

  scale_linetype_manual(values = c(
    "Diagonalisation" = "dashed",
    "WKB" = "dotted"
  )) +

  labs(
    title = "Splitting Energy vs 1/(a^2)",
    linetype = "Theory"
  )

ggplot(EA_data, aes(x = 1/(a * a), y = E0MCMC, color = factor(Separation))) +
  geom_point() +
  geom_line() +
  labs(title = "E0 vs 1/(a * a)")

# EA_data <- EA_data %>%
#   mutate(
#     E0_error = abs(E0MCMC - E0Real),
#     splitting_error = abs((E1 - E0) - (E1Real - E0Real))
#   )

### Most exact wave function vs Potential - Run after running the rest of the code for DWP, then QHO

# Potentials

{
  potential_scale <- 5

  V_QHO <- function(x) { omegaQHO * omegaQHO * mQHO * 0.5 * x * x }
    
  QHO_potential_df <- data.frame(
  x = x_values,
  V_QHO = V_QHO(x_values))

  V_DWP <- function(x) { lambdaDWP * (x * x - wellCentres * wellCentres) * (x * x - wellCentres * wellCentres) }
  
  DWP_potential_df <- data.frame(
  x = x_values,
  V_DWP = V_DWP(x_values))
}

{
  ggplot(wave_df, aes(x = x)) +
    # Exact QHO wave function
    geom_line(aes(y = psiAnalytical, color = "QHO Wave function"), linetype = "dashed") +

    # Diagonalisation DWP wave function
    geom_line(data = diagWFData, aes(x = x, y = psi0, color = "DWP Wave function"), linetype = "dashed") +

    # QHO Potential
    geom_line(data = QHO_potential_df, aes(x = x_values, y = V_QHO / potential_scale, color = "QHO Potential")) +

    # Double Well Potential
    geom_line(data = DWP_potential_df, aes(x = x_values, y = V_DWP / potential_scale, color = "DWP Potential")) +

    labs(title = "QHO and DWP Wave Functions and Potentials", x = "Position x", y = "Psi") +

    coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 0.8)) + 

    scale_y_continuous(
      name = expression("Wavefunction " ~ psi[0](x)),
      sec.axis = sec_axis(~ . * potential_scale,  name = "Potential V(x)")) +

    scale_color_manual(
      name = "",
      values = c("DWP Potential" = "black", "DWP Wave function" = "black", "QHO Potential" = "red", "QHO Wave function" = "red")) + 

    theme_minimal(base_size = 14)
  ggsave("Figures/Wave_functions_and_potentials.png", width = 8, height = 4, dpi = 600)
}

# nolint end