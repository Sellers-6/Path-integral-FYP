# nolint start

{
  library(ggplot2)
  library(gridExtra)  # side by side plots
  library(dplyr)
  library(rhdf5)
  library(minpack.lm) # better exponential fitting
  library(expm) # matrix exponentials
}

#  System
{
  # sys <- "QHO"
  sys <- "DWP"
}

# Overwrite figures?
save_png <- FALSE

# Read data
{
  dataFile <- "dataTESTING.h5"
  
  E0ThermData <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/E0Therm"))))
  accRateThermData <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/accRateTherm"))))
  E0Data <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/E0"))))
  accRateData <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/accRate"))))
  histogramData <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/histogram"))))
  instantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/instantons"))))
  antiInstantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/antiInstantons"))))
  Gx1x1Data <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/Gx1x1"))))
  Gx1x2Data <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/Gx1x2"))))
  Gx1x3Data <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/Gx1x3"))))
  Gx2x2Data <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/Gx2x2"))))
  Gx2x3Data <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/Gx2x3"))))
  Gx3x3Data <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/Gx3x3"))))
  headerData <- as.numeric(unlist(h5read(dataFile, paste0("/", sys, "/headerInfo"))))
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
    if (save_png) {
      ggsave("Figures/QHO_thermalisation.png", width = 6, height = 4, dpi = 1200)
    }
  } else if (sys == "DWP") {
    ggplot(data.frame(sweep = sweepIndex, E0 = E0ThermAvgs), aes(x = sweep * thermInterval, y = E0)) +
      geom_line(color = "blue", linewidth = 0.6) +
      labs(x = "Sweep", y = expression("Ground State Energy" ~ E[0]),
        title = "Average Thermalisation of the DWP") +
      theme_minimal(base_size = 14)
    if (save_png) {
      ggsave("Figures/DWP_thermalisation.png", width = 6, height = 4, dpi = 1200)
    }
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
      geom_line(color = "blue", linewidth = 0.6) +
      labs(x = "Measure", y = expression("Ground State Energy" ~ E[0]), 
        title = "Average Decorrelation of the QHO") +
      theme_minimal(base_size = 14)
    if (save_png) {
      ggsave("Figures/QHO_decorrelation.png", width = 6, height = 4, dpi = 1200)  
    }
  } else if (sys == "DWP") {
    ggplot(data.frame(sweep = measureIndex, E0 = E0Data[1:decorrCheck]), aes(x = sweep, y = E0)) + 
      geom_line(color = "blue", linewidth = 0.6) +
      labs(x = "Measure", y = expression("Ground State Energy" ~ E[0]), 
        title = "Average Decorrelation of the QHO") +
      theme_minimal(base_size = 14)
    if (save_png) {
      ggsave("Figures/DWP_decorrelation.png", width = 6, height = 4, dpi = 1200)  
    }
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
  if (save_png) {
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
  else {
    acfResult <- acf(E0Data[1:measures], lag.max = measures - 1)
  }
}

# Calculating E0 and its error
{
  E0Split <- split(E0Data, rep(1:repeats, each = measures)) # Split E0 into repeats
  
  E0RepeatAvg <- sapply(E0Split, mean)  # Take mean per repeat

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

# Histogram and shapiro test
{
  bins <- 7   # Can be increased for high number of repeats

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
    if (save_png) {
      png("Figures/QHO_Histogram_qq.png", width = 1200, height = 800)
          grid.arrange(histPlot, qqPlot, ncol = 2)
          dev.off()
    }
    else {
      grid.arrange(histPlot, qqPlot, ncol = 2)
    }
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
    if (save_png) {
      png("Figures/DWP_Histogram_qq.png", width = 1200, height = 800)
          grid.arrange(histPlot, qqPlot, ncol = 2)
          dev.off()
    }
    else {
      grid.arrange(histPlot, qqPlot, ncol = 2)
    }
  }
}

# histPlot # Show histogram

# qqPlot # Show QQ plot

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
    if (save_png) {
      ggsave("Figures/QHO_wave_functions.png", width = 8, height = 4, dpi = 1200)  
    }
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
    if (save_png) {
      ggsave("Figures/DWP_wave_functions.png", width = 8, height = 4, dpi = 1200)
    }
  }
}

### x1 x1 correlation function

# Exponential fit

{
  noiseless_region_x1x1 <- 50 # Adjust to the length of the noiseless region
  dfCorr <- data.frame(
    lag = 0:(noiseless_region_x1x1 - 1),
    correlation = Gx1x1Data[1:noiseless_region_x1x1]
  )
  fit_data <- subset(dfCorr, lag >= 0 & lag <= noiseless_region_x1x1)
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
      labs(title = "QHO x1 x1 Correlator with Exponential fit", y = "Two-Point Correlator G(t)", x = "Lag t") +
      scale_color_manual(
        name = "",
        values = c("MCMC Correlator" = "black", "Exponential fit" = "red"))
    if (save_png) {
      ggsave("Figures/QHO_x1x1_correlator.png", width = 7, height = 4, dpi = 1200)
    }
  } else if (sys == "DWP") {
    ggplot(dfCorr, aes(x = lag * latticeSpacing)) +
      geom_point(aes(y = correlation, color = "MCMC Correlator"), size = 0.5) +
      geom_line(aes(y = fit, color = "Exponential fit")) +
      labs(title = "DWP x1 x1 Correlator with Exponential fit", y = "Two-Point Correlator G(t)", x = "Lag t") +
      scale_color_manual(
        name = "",
        values = c("MCMC Correlator" = "black", "Exponential fit" = "red"))
    if (save_png) {
      ggsave("Figures/DWP_x1x1_correlator.png", width = 7, height = 4, dpi = 1200)
    }
  }
}
ggplot(dfCorr, aes(x = lag * latticeSpacing)) +
      geom_point(aes(y = correlation, color = "MCMC Correlator"), size = 0.5) +
      geom_line(aes(y = fit, color = "Exponential fit")) +
      labs(title = "DWP x1 x1 Correlator with Exponential fit", y = "Two-Point Correlator G(t)", x = "Lag t") +
      scale_color_manual(
        name = "",
        values = c("MCMC Correlator" = "black", "Exponential fit" = "red"))

# x2 x2 correlation function
{
  noiseless_region_x2x2 <- 50 # Adjust to the length of the noiseless region
  dfCorr <- data.frame(
    lag = 0:(noiseless_region_x2x2 - 1),
    correlation = Gx2x2Data[0:noiseless_region_x2x2]
  )
}

ggplot(dfCorr, aes(x = lag, y = correlation)) +
  geom_line(color = "#000000") +
  labs(
    title = paste("x2 x2 correlation function for", sys),
    x = "Time (index of the path)",
    y = "G_4(t, 0)"
  )

### Excited energy states

# E1

{
  successfulCounts <- 0; E1 <- 0

  LHS <- 1; RHS <- noiseless_region_x1x1

  correlatorRatios <- numeric(RHS - LHS + 1) # To store the log ratios for each lag

  for (i in LHS:RHS) {
    if (Gx1x1Data[i] <= 0 || Gx1x1Data[i + 1] <= 0) {
      message("Correlation function has non-positive values, cannot compute log ratio.")
    } 
    else {
      E1 <- E1 + E0 + log(Gx1x1Data[i] / Gx1x1Data[i + 1]) / latticeSpacing
      correlatorRatios[i - LHS + 1] <- log(Gx1x1Data[i] / Gx1x1Data[i + 1]) / latticeSpacing
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
    if (save_png) {
      ggsave("Figures/QHO_log_ratio.png", width = 6, height = 4, dpi = 1200)
    }

  } else if (sys == "DWP") {
    ggplot(data.frame(lag = 1:(length(correlatorRatios)), ratio = correlatorRatios), aes(x = lag * latticeSpacing, y = ratio)) +
      geom_line(aes(color = "Log Ratio")) +
      geom_hline(data = data.frame(y = E1 - E0), aes(yintercept = y, color = "Splitting Energy"), linetype = "dashed") +
      labs(title = "Finding the splitting energy for the DWP",
          x = "Lag t", y = expression("Log(" ~ G[2](t) ~ ") / Log(" ~ G[2](t+a) ~ ")")) +
      scale_color_manual(name = "", values = c("Log Ratio" = "blue", "Splitting Energy" = "red"))
    if (save_png) {
      ggsave("Figures/DWP_log_ratio.png", width = 6, height = 4, dpi = 1200)
    }
  } 
}

ggplot(data.frame(lag = 1:(length(correlatorRatios)), ratio = correlatorRatios), aes(x = lag * latticeSpacing, y = ratio)) +
      geom_line(aes(color = "Log Ratio")) +
      geom_hline(data = data.frame(y = E1 - E0), aes(yintercept = y, color = "Splitting Energy"), linetype = "dashed") +
      labs(title = "Finding the splitting energy for the QHO",
          x = "Lag t", y = expression("Log(" ~ G[2](t) ~ ") / Log(" ~ G[2](t+a) ~ ")")) +
      scale_color_manual(name = "", values = c("Log Ratio" = "blue", "Splitting Energy" = "red"))

E1
E1 - E0
Split_diag

{
  successfulCounts <- 0; E2 <- 0

  for (i in 1:noiseless_region_x2x2) {
    if (Gx2x2Data[i] <= 0 || Gx2x2Data[i + 1] <= 0) {
      message("Correlation function has non-positive values, cannot compute log ratio.")
    } 
    else {
      E2 <- E2 + mean(E0RepeatAvg) + log(Gx2x2Data[i] / Gx2x2Data[i + 1]) / latticeSpacing
      successfulCounts <- successfulCounts + 1
    }
  }
  E2 <- E2 / successfulCounts; E2
}

E2 - E0

### GEVP - Alternative method of extracting excited energy state.division

{
  # Matrix of correlators with vector entries "pathLength" long
  C <- array(0, dim = c(3, 3, pathLength)) 

  C[1,1,] <- Gx1x1Data
  C[1,2,] <- Gx1x2Data
  C[1,3,] <- Gx1x3Data

  C[2,1,] <- Gx1x2Data
  C[2,2,] <- Gx2x2Data
  C[2,3,] <- Gx2x3Data

  C[3,1,] <- Gx1x3Data
  C[3,2,] <- Gx2x3Data
  C[3,3,] <- Gx3x3Data

  t0 <- 1

  eigvals <- matrix(0, nrow = pathLength, ncol = 3)

  C0 <- C[,,t0]
  C0 <- (C0 + t(C0)) / 2

  C0_inv_sqrt <- solve(sqrtm(C0))

  for (t in 1:pathLength) {
    Ct <- C[,,t]
    Ct <- (Ct + t(Ct)) / 2
    
    M <- C0_inv_sqrt %*% Ct %*% C0_inv_sqrt
    gevp <- eigen(M)
    
    eigvals[t, ] <- Re(gevp$values)
  }

  eigvals <- t(apply(eigvals, 1, sort, decreasing = TRUE))

  # ---- Effective energies from GEVP eigenvalues ----
  E <- matrix(NA, nrow = pathLength - 1, ncol = 3)

  for (t in 1:(pathLength - 1)) {

    ratio <- eigvals[t, ] / eigvals[t + 1, ]

    # avoid invalid values
    ratio[ratio <= 0] <- NA

    E[t, ] <- -log(ratio) / latticeSpacing
  }
}

{
  # Matrix of correlators with vector entries "pathLength" long
  C <- array(0, dim = c(2, 2, pathLength)) 

  C[1,1,] <- Gx1x1Data
  C[1,2,] <- Gx1x2Data

  C[2,1,] <- Gx1x2Data
  C[2,2,] <- Gx2x2Data

  t0 <- 2

  eigvals <- matrix(0, nrow = pathLength, ncol = 2)

  C0 <- C[,,t0]
  C0 <- (C0 + t(C0)) / 2

  C0_inv_sqrt <- solve(sqrtm(C0))

  for (t in 1:pathLength) {
    Ct <- C[,,t]
    Ct <- (Ct + t(Ct)) / 2
    
    M <- C0_inv_sqrt %*% Ct %*% C0_inv_sqrt
    gevp <- eigen(M)
    
    eigvals[t, ] <- Re(gevp$values)
  }

  eigvals <- t(apply(eigvals, 1, sort, decreasing = TRUE))

  # ---- Effective energies from GEVP eigenvalues ----
  E <- matrix(NA, nrow = pathLength - 1, ncol = 2)

  for (t in 1:(pathLength - 1)) {

    ratio <- eigvals[t, ] / eigvals[t + 1, ]

    # avoid invalid values
    ratio[ratio <= 0] <- NA

    E[t, ] <- -log(ratio) / latticeSpacing
  }
}

find_flat_region <- function(x, window = 10) {
  
  n <- length(x)
  
  if (window > n) stop("window is larger than data length")
  
  best_var <- Inf
  best_start <- 1
  
  # slide window across data
  for (i in 1:(n - window + 1)) {
    
    segment <- x[i:(i + window - 1)]
    
    # ignore NA values safely
    segment <- segment[!is.na(segment)]
    
    if (length(segment) < 2) next
    
    v <- var(segment)
    
    if (v < best_var) {
      best_var <- v
      best_start <- i
    }
  }
  
  best_region <- x[best_start:(best_start + window - 1)]
  
  list(
    mean = mean(best_region, na.rm = TRUE),
    start_index = best_start,
    end_index = best_start + window - 1,
    variance = best_var,
    region = best_region
  )
}

cut_off <- 50

df <- data.frame(
  t = 1:cut_off,
  E0 = E[,1][1:cut_off]
)

plateau_E0 <- find_flat_region(E[,1], window = 10)
plateau_E0$mean

ggplot(df, aes(x = t, y = E0)) +
  geom_line() +
  labs(x = "t", y = "E0(t)", title = "Effective Energy (Ground state)")

df <- data.frame(
  t = 1:cut_off,
  E1 = E[,2][1:cut_off]
)

plateau_E1 <- find_flat_region(E[,2], window = 10)
plateau_E1$mean

ggplot(df, aes(x = t, y = E1)) +
  geom_line() +
  labs(x = "t", y = "E1(t)", title = "Effective Energy (1st excited state)")

df <- data.frame(
  t = 1:cut_off,
  E2 = E[,3][1:cut_off]
)

plateau_E2 <- find_flat_region(E[,3], window = 10)
plateau_E2$mean

ggplot(df, aes(x = t, y = E2)) +
  geom_line() +
  labs(x = "t", y = "E2(t)", title = "Effective Energy (2nd excited state)")



Split_diag

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

#####################################################################################################################################
####################################################### Read new data ###############################################################
#####################################################################################################################################

{
  library(ggplot2)
  library(gridExtra)  # side by side plots
  library(dplyr)
  library(rhdf5)
  library(minpack.lm) # better exponential fitting
  library(expm) # eigenvalue solving
}

{
  dataFile <- "data.h5"
  sys <- "DWP"

  well_centres_vec <- c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
  lattice_spacings_vec <- c(0.1)
  epsilon_vec <- c(0.4)
  path_lengths <- c(100000)
  betas <- c(10000)

  wellReps <- length(well_centres_vec)
  latReps <- length(lattice_spacings_vec)

  for (i in seq_along(lattice_spacings_vec)) {
    if (lattice_spacings_vec[i] * path_lengths[i] != betas[i]) {
      stop("Check beta vector")
    }
  }

  # Build paths
  full_paths <- matrix("", nrow = wellReps, ncol = latReps)

  for (i in seq_along(well_centres_vec)) {
    wc_str <- sprintf("%.6f", well_centres_vec[i])
    
    for (j in seq_along(lattice_spacings_vec)) {
      ls_str <- sprintf("%.6f", lattice_spacings_vec[j])
      
      full_paths[i, j] <- paste0("/", sys, "/", wc_str, "/", ls_str, "/")
    }
  }

  # Read header data
  headerSize <- 12

  headerData <- array(NA_real_, dim = c(wellReps, latReps, headerSize))

  for (i in seq_along(well_centres_vec)) {
    for (j in seq_along(lattice_spacings_vec)) {
      path <- paste0(full_paths[i, j], "headerInfo")
      
      tmp <- as.numeric(unlist(h5read(dataFile, path)))
      
      if (length(tmp) != headerSize) {
        stop(paste("Header size mismatch at", path))
      }
      
      headerData[i, j, ] <- tmp
    }
  }

  # Extract constants (these are same for all runs)
  decorrSweeps    <- headerData[1, 1, 5]
  measures        <- headerData[1, 1, 8]
  repeats         <- headerData[1, 1, 9]
  numBins         <- headerData[1, 1, 10]
  lambdaDWP       <- headerData[1, 1, 12]


  totalMeasures <- measures * repeats
  totalAcc <- measures * repeats

  # Allocate data arrays
  E0Data             <- array(0, dim = c(wellReps, latReps, totalMeasures))
  accRateData        <- array(0, dim = c(wellReps, latReps, totalAcc))
  histogramData      <- array(0, dim = c(wellReps, latReps, numBins))
  instantonsData     <- array(0, dim = c(wellReps, latReps, totalMeasures))
  antiInstantonsData <- array(0, dim = c(wellReps, latReps, totalMeasures))

  # Correlators
  corrSize <- path_lengths[1]

  Gx1x1Data <- array(0, dim = c(wellReps, latReps, corrSize))
  Gx1x2Data <- array(0, dim = c(wellReps, latReps, corrSize))
  Gx1x3Data <- array(0, dim = c(wellReps, latReps, corrSize))
  Gx2x2Data <- array(0, dim = c(wellReps, latReps, corrSize))
  Gx2x3Data <- array(0, dim = c(wellReps, latReps, corrSize))
  Gx3x3Data <- array(0, dim = c(wellReps, latReps, corrSize))

  # Read data
  for (i in seq_along(well_centres_vec)) {
    for (j in seq_along(lattice_spacings_vec)) {
      
      base <- full_paths[i, j]
      
      E0Data[i, j, ]             <- as.numeric(unlist(h5read(dataFile, paste0(base, "E0"))))
      accRateData[i, j, ]        <- as.numeric(unlist(h5read(dataFile, paste0(base, "accRate"))))
      histogramData[i, j, ]      <- as.numeric(unlist(h5read(dataFile, paste0(base, "histogram"))))
      instantonsData[i, j, ]     <- as.numeric(unlist(h5read(dataFile, paste0(base, "instantons"))))
      antiInstantonsData[i, j, ] <- as.numeric(unlist(h5read(dataFile, paste0(base, "antiInstantons"))))
      
      Gx1x1Data[i, j, ] <- as.numeric(unlist(h5read(dataFile, paste0(base, "Gx1x1"))))
      Gx1x2Data[i, j, ] <- as.numeric(unlist(h5read(dataFile, paste0(base, "Gx1x2"))))
      Gx1x3Data[i, j, ] <- as.numeric(unlist(h5read(dataFile, paste0(base, "Gx1x3"))))
      Gx2x2Data[i, j, ] <- as.numeric(unlist(h5read(dataFile, paste0(base, "Gx2x2"))))
      Gx2x3Data[i, j, ] <- as.numeric(unlist(h5read(dataFile, paste0(base, "Gx2x3"))))
      Gx3x3Data[i, j, ] <- as.numeric(unlist(h5read(dataFile, paste0(base, "Gx3x3"))))
    }
  }
}

{
  # Theoretical calculations (vectorised over well separations)
  wellCentres <- well_centres_vec

  omegaDWP <- sqrt(8 * lambdaDWP * wellCentres^2)
  S_inst   <- (2/3) * omegaDWP * (wellCentres^2)

  alpha <- 1 / 12
  K <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5)

  Split_WKB <- K * exp(-S_inst)
  E0_WKB <- 0.5 * omegaDWP - (Split_WKB / 2)
  E1_WKB <- 0.5 * omegaDWP + (Split_WKB / 2)

  # Diagonalisation data
  diagEnergiesData <- read.csv("DWP diagonalisation/DWP diagonalisation/energies.csv")

  energy_row <- as.integer(round((wellCentres - 0.9) * 10))

  E0_diag    <- diagEnergiesData[energy_row, 2]
  E1_diag    <- diagEnergiesData[energy_row, 3]
  Split_diag <- diagEnergiesData[energy_row, 4]
}

# Check acceptance rate of any iteration by changing the indices of accRateData
mean(accRateData[11, 1, ]) * 100

E0_vals <- matrix(0, nrow = wellReps, ncol = latReps)
E0_errors  <- matrix(0, nrow = wellReps, ncol = latReps)

for (i in seq_len(wellReps)) {
  for (j in seq_len(latReps)) {
    
    # Extract time series for this parameter set
    E0_series <- E0Data[i, j, ]
    
    # Split into repeats
    E0Split <- split(E0_series, rep(1:repeats, each = measures))
    
    # Mean per repeat
    E0RepeatAvg <- sapply(E0Split, mean)
    
    # Final estimates
    E0_vals[i, j] <- mean(E0RepeatAvg)
    E0_errors[i, j]  <- sd(E0RepeatAvg) / sqrt(length(E0RepeatAvg))
  }
}

E0_vals[,1]
E0_errors[,1]

E0_vals[,1]; E0_vals[,1] + E0_errors[,1]; E0_vals[,1] - E0_errors[,1]

# Correlators - Naive version

fit_correlator <- function(corr, latticeSpacing, noiseless_region, DeltaE_guess) {
  
  lag <- 0:(noiseless_region - 1)
  corr_slice <- corr[1:noiseless_region]
  
  df <- data.frame(lag = lag, correlation = corr_slice)
  
  fit <- nlsLM(
    correlation ~ A * exp(-DeltaE * lag * latticeSpacing) + c,
    data = df,
    start = list(
      A = corr_slice[1],
      c = min(corr_slice),
      DeltaE = DeltaE_guess
    ),
    control = nls.lm.control(maxiter = 1024)
  )
  
  list(
    DeltaE = coef(fit)["DeltaE"],
    fit = predict(fit, newdata = df),
    df = df
  )
}

log_ratio_energy <- function(corr, latticeSpacing, E0, max_lag) {
  
  ratios <- numeric(max_lag)
  valid <- 0
  
  for (i in 1:max_lag) {
    if (corr[i] > 0 && corr[i + 1] > 0) {
      ratios[i] <- log(corr[i] / corr[i + 1]) / latticeSpacing
      valid <- valid + 1
    } else {
      ratios[i] <- NA
    }
  }
  
  ratios <- ratios[!is.na(ratios)]
  
  DeltaE <- mean(ratios)
  E <- E0 + DeltaE
  
  list(
    E = E,
    DeltaE = DeltaE,
    ratios = ratios
  )
}

E1_vals <- matrix(NA, wellReps, latReps)
E1_err  <- matrix(NA, wellReps, latReps)

E2_vals <- matrix(NA, wellReps, latReps)

noiseless_region <- 10

for (i in seq_len(wellReps)) {
  for (j in seq_len(latReps)) {
    
    latticeSpacing <- lattice_spacings_vec[j]
    E0_here <- E0_vals[i, j]
    
    # -------------------------
    # x1x1 → E1
    # -------------------------
    corr1 <- Gx1x1Data[i, j, ]
    
    res1 <- log_ratio_energy(
      corr = corr1,
      latticeSpacing = latticeSpacing,
      E0 = E0_here,
      max_lag = noiseless_region
    )
    
    E1_vals[i, j] <- res1$E
    
    # crude error estimate
    E1_err[i, j] <- sd(res1$ratios) / sqrt(length(res1$ratios))
    
    
    # -------------------------
    # x2x2 → E2
    # -------------------------
    corr2 <- Gx2x2Data[i, j, ]
    
    res2 <- log_ratio_energy(
      corr = corr2,
      latticeSpacing = latticeSpacing,
      E0 = E0_here,
      max_lag = noiseless_region
    )
    
    E2_vals[i, j] <- res2$E
  }
}

Split_fit <- matrix(NA, wellReps, latReps)

for (i in seq_len(wellReps)) {
  for (j in seq_len(latReps)) {
    
    latticeSpacing <- lattice_spacings_vec[j]
    corr <- Gx1x1Data[i, j, ]
    
    fit_res <- fit_correlator(
      corr,
      latticeSpacing,
      noiseless_region,
      DeltaE_guess = 0.1 # Put the diagonalisation splitting energy here as a guess
    )
    
    Split_fit[i, j] <- fit_res$DeltaE
  }
}

plot_correlator <- function(corr, latticeSpacing, fit_res = NULL) {
  
  n <- length(corr)
  df <- data.frame(
    lag = 0:(n - 1),
    corr = corr
  )
  
  p <- ggplot(df, aes(x = lag * latticeSpacing, y = corr)) +
    geom_line(size = 0.5) +
    labs(x = "t", y = "G(t)")
  
  if (!is.null(fit_res)) {
    df$fit <- fit_res$fit
    p <- p + geom_line(aes(y = fit), color = "red")
  }
  
  p
}

plot_correlator(Gx1x1Data[10,1,][1:10], 0.1)

E1_vals[,1]

E2_vals[,1]

Split_fit[,1]

#################################################
# GEVP method for finding excited energy states #
#################################################

build_C_matrix <- function(Glist, pathLength) {
  
  if (length(Glist) == 3) {
    n <- 2
  } else if (length(Glist) == 6) {
    n <- 3
  } else {
    stop("Unsupported basis size")
  }
  
  C <- array(0, dim = c(n, n, pathLength))
  
  if (n == 2) {
    C[1,1,] <- Glist[[1]]
    C[1,2,] <- Glist[[2]]
    C[2,1,] <- Glist[[2]]
    C[2,2,] <- Glist[[3]]
    
  } else {
    C[1,1,] <- Glist[[1]]
    C[1,2,] <- Glist[[2]]
    C[1,3,] <- Glist[[3]]
    
    C[2,1,] <- Glist[[2]]
    C[2,2,] <- Glist[[4]]
    C[2,3,] <- Glist[[5]]
    
    C[3,1,] <- Glist[[3]]
    C[3,2,] <- Glist[[5]]
    C[3,3,] <- Glist[[6]]
  }
  
  return(C)
}

solve_gevp <- function(C, latticeSpacing, t0 = 1) {
  
  pathLength <- dim(C)[3]
  n <- dim(C)[1]
  
  eigvals <- matrix(NA, nrow = pathLength, ncol = n)
  
  # Reference matrix
  C0 <- (C[,,t0] + t(C[,,t0])) / 2
  C0_inv_sqrt <- solve(sqrtm(C0))
  
  for (t in 1:pathLength) {
    
    Ct <- (C[,,t] + t(C[,,t])) / 2
    
    M <- C0_inv_sqrt %*% Ct %*% C0_inv_sqrt
    
    ev <- eigen(M, symmetric = TRUE)
    eigvals[t, ] <- sort(Re(ev$values), decreasing = TRUE)
  }
  
  # Effective energies
  E <- matrix(NA, nrow = pathLength - 1, ncol = n)
  
  for (t in 1:(pathLength - 1)) {
    ratio <- eigvals[t + 1, ] / eigvals[t, ]
    ratio[ratio <= 0] <- NA
    E[t, ] <- -log(ratio) / latticeSpacing
  }
  
  E
}

find_plateau <- function(x, window = 10) {
  n <- length(x)
  best_var <- Inf
  best_start <- 1
  
  for (i in 1:(n - window + 1)) {
    seg <- x[i:(i + window - 1)]
    seg <- seg[!is.na(seg)]
    if (length(seg) < 2) next
    
    v <- var(seg)
    if (v < best_var) {
      best_var <- v
      best_start <- i
    }
  }
  
  region <- x[best_start:(best_start + window - 1)]
  
  list(
    mean = mean(region, na.rm = TRUE),
    start = best_start,
    end = best_start + window - 1,
    var = best_var
  )
}

plot_effective_energy <- function(E, latticeSpacing, state = 1, cutoff = 25) {
  
  vals <- E[, state]
  
  # remove NA values
  valid_idx <- which(!is.na(vals))
  
  if (length(valid_idx) < 2) {
    warning("Not enough valid points to plot")
    return(NULL)
  }
  
  # restrict to cutoff safely
  idx <- valid_idx[valid_idx <= cutoff]
  
  df <- data.frame(
    t = idx * latticeSpacing,
    E = vals[idx]
  )
  
  ggplot(df, aes(x = t, y = E)) +
    geom_line() +
    labs(
      x = "τ",
      y = paste0("E", state - 1, "(τ)"),
      title = paste("Effective Energy (state", state - 1, ")")
    )
}

E0_gevp <- matrix(NA, wellReps, latReps)
E1_gevp <- matrix(NA, wellReps, latReps)
E2_gevp <- matrix(NA, wellReps, latReps)

pathLength <- dim(Gx1x1Data)[3]   # safer than recomputing
nStates <- 3

E_array <- array(NA_real_,
                 dim = c(wellReps, latReps, pathLength - 1, nStates))

for (i in seq_len(wellReps)) {
  for (j in seq_len(latReps)) {
    
    latticeSpacing <- lattice_spacings_vec[j]
    
    Glist <- list(
      Gx1x1Data[i,j,],
      Gx1x2Data[i,j,],
      Gx1x3Data[i,j,],
      Gx2x2Data[i,j,],
      Gx2x3Data[i,j,],
      Gx3x3Data[i,j,]
    )
    
    C <- build_C_matrix(Glist, pathLength)
    E <- solve_gevp(C, latticeSpacing, t0 = 2)
    
    # ✅ store full time-dependent energies
    E_array[i, j, , ] <- E
    
    # plateau extraction (unchanged)
    p0 <- find_plateau(E[,1])
    p1 <- find_plateau(E[,2])
    p2 <- find_plateau(E[,3])
    
    E0_gevp[i,j] <- p0$mean
    E1_gevp[i,j] <- p1$mean
    E2_gevp[i,j] <- p2$mean
  }
}

i <- 10
j <- 1
latticeSpacing <- lattice_spacings_vec[j]

E <- E_array[i, j, , ]

plot_effective_energy(E, latticeSpacing, state = 1)
plot_effective_energy(E, latticeSpacing, state = 2)
plot_effective_energy(E, latticeSpacing, state = 3)

E0_gevp[,1]

E <- solve_gevp(C, latticeSpacing, t0 = 2)

plot_effective_energy(E, latticeSpacing, state = 1)

E0_gevp

# Create graph of splitting energy vs well separation with error bars

diagEnergiesData <- read.csv("DWP diagonalisation/DWP diagonalisation/energies.csv")
diagWFData <- read.csv("DWP diagonalisation/DWP diagonalisation/wavefunctions.csv")

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
  geom_hline(data = theory_df,
             aes(yintercept = E0WKB, linetype = "WKB"),
             color = "black") +

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

# nolint end