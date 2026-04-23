# nolint start

#####################################################################################################################################
################################################ Quantum Harmonic Oscillator ########################################################
#####################################################################################################################################

###~~~~~~~~~~~###
### Libraries ###
###~~~~~~~~~~~###

{
  library(ggplot2)
  library(gridExtra)    # Side by side plots
  library(dplyr)
  library(rhdf5)        # Read h5 files
  library(minpack.lm)   # Better exponential fitting
  library(expm)         # Matrix exponentials
}

# Overwrite figures?
save_png <- FALSE

###~~~~~~~~~~~###
### Read data ###
###~~~~~~~~~~~###

{
  dataFile <- "data.h5"
  
  accRateThermData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/accRateTherm"))))
  accRateData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/accRate"))))
  E0ThermData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/E0Therm"))))
  E0Data <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/E0"))))
  histogramData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/histogram"))))
  Gx1x1Data <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/Gx1x1"))))
  Gx2x2Data <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/Gx2x2"))))
  instantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/instantons"))))
  antiInstantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/antiInstantons"))))
  headerData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/headerInfo"))))
}

###~~~~~~~~~~~~~~~~~~~~###
### Variables from C++ ###
###~~~~~~~~~~~~~~~~~~~~###

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
  mQHO <- headerData[11]
  omegaQHO <- headerData[12]
  
  beta <- pathLength * latticeSpacing
  thermMeasures <- thermSweeps / thermInterval
}

###~~~~~~~~~~~~~~~~~###
### Acceptance rate ###
###~~~~~~~~~~~~~~~~~###

mean(accRateThermData) * 100 # Should be between 50 and 80%

###~~~~~~~~~~~~~~~~###
### Thermalisation ###
###~~~~~~~~~~~~~~~~###

{
  # Split into repeats, each column is one repeat
  E0Mat <- matrix(E0ThermData, nrow = thermMeasures, ncol = repeats, byrow = FALSE)

  # Compute average thermalisation across repeats
  E0ThermAvgs <- rowMeans(E0Mat)

  # Sweep index for plotting
  sweepIndex <- seq_len(thermMeasures)

  # Data frame for plotting
  E0ThermDF <- data.frame(sweep = sweepIndex, E0 = E0ThermAvgs)
  
  # Plot
  ggplot(data.frame(sweep = sweepIndex, E0 = E0ThermAvgs), aes(x = sweep * thermInterval, y = E0)) +
      geom_line(color = "blue", linewidth = 0.6) +
      labs(x = "Sweep", y = expression("Ground State Energy" ~ E[0]), title = "Average Thermalisation of the QHO") +
      theme_minimal(base_size = 14)
}

if (save_png) { ggsave("Figures/QHO_thermalisation.png", width = 6, height = 4, dpi = 1200) }

###~~~~~~~~~~~~~~~###
### Decorrelation ###
###~~~~~~~~~~~~~~~###

{
  if (measures > 50) {
    decorrCheck <- 50
  } else {
    decorrCheck <- measures
  }
  measureIndex <- seq_len(decorrCheck)
  ggplot(data.frame(sweep = measureIndex, E0 = E0Data[1:decorrCheck]), aes(x = sweep, y = E0)) + 
    geom_line(color = "blue", linewidth = 0.6) +
    labs(x = "Measure", y = expression("Ground State Energy" ~ E[0]), title = "Average Decorrelation of the QHO") +
    theme_minimal(base_size = 14)
}

if (save_png) { ggsave("Figures/QHO_decorrelation.png", width = 6, height = 4, dpi = 1200) }

###~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Auto-correlation function ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~###

{
  acfResult <- acf(E0Data[1:measures], lag.max = measures - 1, plot = TRUE)
  acfVals <- acfResult$acf
  mean(abs(as.vector(acfVals)[2:length(acfVals)])) # For decorrelated data, we want low acf values (roughly < 0.1)
}

if (save_png) {
    png("Figures/acf_plot_QHO.png", width = 800, height = 600)

    acfResult <- acf(E0Data[1:measures], lag.max = measures - 1, plot = FALSE)

    plot(acfResult, main = "Autocorrelation of Ground State Energy for the QHO", xlab = "Lag (sweeps)", ylab = "Autocorrelation")

    dev.off()
}

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Calculate ground state energy ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

{
  E0Split <- split(E0Data, rep(1:repeats, each = measures)) # Split E0 into repeats
  
  E0RepeatAvg <- sapply(E0Split, mean)  # Take mean per repeat

  E0 <- mean(E0RepeatAvg)

  E0StandardError <- sd(E0RepeatAvg) / sqrt(length(E0RepeatAvg))

  E0_pc_error <- abs(((E0 - 0.5) / 0.5) * 100) # Percent error
}

E0; mean(E0RepeatAvg) + E0StandardError; mean(E0RepeatAvg) - E0StandardError ; E0_pc_error; E0StandardError

###~~~~~~~~~~~~~~~~~###
### Normality tests ###
###~~~~~~~~~~~~~~~~~###

{
  bins <- 7

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
    geom_line(aes(x = continuousE0, y = normDist), color = "red", linewidth = 1) +
    labs(title = "Ground State Energy Histogram for the QHO",
        x = expression("Ground State Energy " ~ E[0]), y = "Probability Density") 

  shapiroTest <- shapiro.test(E0RepeatAvg) # Run a Shapiro test on the data, values of p < 0.05 are not normally distributed
  qqPlot <- ggplot(data.frame(E0RepeatAvg), aes(sample = E0RepeatAvg)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = paste0("QQ Plot (Shapiro-Wilk p = ", round(shapiroTest$p.value, 4), ")"),
        x = "Theoretical Quantiles", y = "Sample Quantiles")

  grid.arrange(histPlot, qqPlot, ncol = 2)
}

if (save_png) {
      png("Figures/QHO_Histogram_qq.png", width = 1200, height = 800)
          grid.arrange(histPlot, qqPlot, ncol = 2)
          dev.off()
    }

# histPlot # Show histogram

# qqPlot # Show QQ plot

###~~~~~~~~~~~~~~~###
### Wave function ###
###~~~~~~~~~~~~~~~###

# Histogram data frame creation

{
  # Histogram range
  sigmaQHO <- 1 / sqrt(2 * mQHO * omegaQHO)
  xMax <- ceiling(4 * sigmaQHO) + 1  # 4 standard deviations of analytic ground state plus padding
  xMin <- -xMax

  binWidth <- (xMax - xMin) / numBins
  x_values <- seq(xMin + binWidth/2, xMax - binWidth/2, length.out = numBins)

  hist_matrix <- matrix(histogramData, nrow = repeats, ncol = numBins, byrow = TRUE)
  hist_avg <- colMeans(hist_matrix)

  prob_density <- hist_avg / (sum(hist_avg) * binWidth)

  hist_df <- data.frame(x = x_values, probability = prob_density)

  psi <- sqrt(prob_density)

  wave_df <- data.frame(x = x_values, psi = psi)

  # Overlay analytical wave function 
  psiAnalytical <- exp(-(x_values^2) / 2)
  
  # Normalise the analytical wave function 
  psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) 

  wave_df$psiAnalytical <- psiAnalytical
}

# Comparing wave functions

{
  ggplot(wave_df, aes(x = x)) +
    # MCMC wave function
    geom_point(aes(y = psi, color = "MCMC"), size = 1.5) +
    geom_line(aes(y = psi, color = "MCMC"), linewidth = 0.7, linetype = "dotted") +

    # Exact wave function
    geom_line(aes(y = psiAnalytical, color = "Exact")) +

    labs(title = paste("Quantum Harmonic Oscillator Ground State Wave Function"), x = "Position", y = expression("Wavefunction " ~ psi[0](x))) +

    scale_color_manual(name = "", values = c("MCMC" = "red", "Exact" = "black")) + 

    theme_minimal(base_size = 14)
}
   
if (save_png) { ggsave("Figures/QHO_wave_functions.png", width = 8, height = 4, dpi = 1200) }

###~~~~~~~~~~~~~~~~~~~~~~~###
### Excited energy states ###
###~~~~~~~~~~~~~~~~~~~~~~~###

### Log ratio method ### 

find_plateau <- function(ratios, width) {
  n <- length(ratios)
  
  # Check if data is shorter than width
  if (n <= width) return(1:n)
  
  variances <- sapply(1:(n - width + 1), function(i) {
    var(ratios[i:(i + width - 1)])
  })
  
  # Find the start index of the window with the minimum variance
  best_start <- which.min(variances)
  
  return(best_start:(best_start + width - 1))
}

log_ratio_energy <- function(correlator, latticeSpacing, E0, plateau_width) {
  # Calculate all log ratios
  max_lag <- length(correlator) / 2
  raw_ratios <- numeric(max_lag)
  for (t in 1:max_lag) {
    if (correlator[t] > 0 && correlator[t + 1] > 0) {
      raw_ratios[t] <- log(correlator[t] / correlator[t + 1]) / latticeSpacing
    } 
    else { raw_ratios[t] <- NA }
  }
  
  clean_ratios <- raw_ratios[!is.na(raw_ratios)]
  
  # Find the flattest window of the specified width
  plateau_indices <- find_plateau(clean_ratios, width = plateau_width)
  plateau_values  <- clean_ratios[plateau_indices]
  
  DeltaE <- mean(plateau_values)
  E      <- E0 + DeltaE
  
  list(
    E = E,
    DeltaE = DeltaE,
    ratios = clean_ratios,
    plateau = plateau_values,
    plateau_indices = plateau_indices
  )
}

# E1

{
  log_ratio_results_1 <- log_ratio_energy(
    correlator = Gx1x1Data,
    latticeSpacing = latticeSpacing,
    E0 = E0,
    plateau_width = 10
  )
      
  E1 <- log_ratio_results_1$E
  delta_E <- log_ratio_results_1$DeltaE
  # Estimate error by the standard deviation of ratios in the best plateau 
  delta_E_err <- sd(log_ratio_results_1$plateau) / sqrt(length(log_ratio_results_1$plateau))

  noiseless_region_1 <- log_ratio_results_1$plateau_indices

  E1_pc_error <- abs(((E1 - 1.5) / 1.5) * 100)
  delta_E_pc_error <- abs(((delta_E - 1.0) / 1.0) * 100)
}

E1; 1.5; E1 - sqrt(delta_E_err^2 + E0StandardError^2); E1 + sqrt(delta_E_err^2 + E0StandardError^2); E1_pc_error; sqrt(delta_E_err^2 + E0StandardError^2)
delta_E; 1.0; delta_E - delta_E_err; delta_E + delta_E_err; delta_E_pc_error

# E2

{
  log_ratio_results_2 <- log_ratio_energy(
    correlator = Gx2x2Data,
    latticeSpacing = latticeSpacing,
    E0 = E0,
    plateau_width = 10
  )
      
  E2 <- log_ratio_results_2$E
  delta_E_2 <- log_ratio_results_2$DeltaE
  # Estimate error by the standard deviation of ratios in the best plateau 
  delta_E_2_err <- sd(log_ratio_results_2$plateau) / sqrt(length(log_ratio_results_2$plateau))

  noiseless_region_2 <- log_ratio_results_2$plateau_indices

  E2_pc_error <- abs(((E2 - 2.5) / 2.5) * 100)
  delta_E_2_pc_error <- abs(((delta_E_2 - 2.0) / 2.0) * 100)
}

E2; 2.5; E2 - delta_E_2_err - E0StandardError; E2 + delta_E_2_err + E0StandardError; E2_pc_error
delta_E_2; 2.0; delta_E_2 - delta_E_2_err; delta_E_2 + delta_E_2_err; delta_E_2_pc_error

### Exponential fit ###

fit_correlator_exp <- function(correlator, noiseless_region, A_guess, DeltaE_guess, latticeSpacing) {
  lag <- (noiseless_region - 1)
  correlator <- correlator[noiseless_region]
  
  dfCorr <- data.frame(
    lag = lag,
    correlation = correlator
  )
  
  c_guess <- tail(correlator, 1)
  
  fit <- nlsLM(
    correlation ~ A * exp(-DeltaE * lag * latticeSpacing) + c,
    data = dfCorr,
    start = list(
      A = A_guess,
      DeltaE = DeltaE_guess,
      c = c_guess),
    control = nls.lm.control(maxiter = 1024)
  )
  
  dfCorr$fit <- predict(fit, newdata = dfCorr)
  
  p <- ggplot(dfCorr, aes(x = lag * latticeSpacing)) +
      geom_point(aes(y = correlation, color = "MCMC Correlator"), size = 1) +
      geom_line(aes(y = fit, color = "Exponential fit"), linewidth = 0.8) +
      labs(title = "Correlator with Exponential Fit",
        subtitle = bquote(Delta * E == .(round(coef(fit)["DeltaE"], 4))), y = "G(t)", x = "t") +
      scale_color_manual(name = "", values = c("MCMC Correlator" = "black", "Exponential fit" = "red")) +
      theme_minimal()
  
  list(
    coefficients = coef(fit),
    fit_object = fit,
    data = dfCorr,
    plot = p
  )
}

# Naively choose noiseless region

{
  Gx1x1Fit <- fit_correlator_exp(correlator = Gx1x1Data, noiseless_region = 1:50, A_guess = 0.5, DeltaE_guess = 1.0, latticeSpacing = latticeSpacing)
  Gx2x2Fit <- fit_correlator_exp(correlator = Gx2x2Data, noiseless_region = 1:30, A_guess = 0.5, DeltaE_guess = 2.0, latticeSpacing = latticeSpacing)
}

Gx1x1Fit$coefficients
print(Gx1x1Fit$plot) 

if (save_png) { ggsave("Figures/QHO_x1x1_correlator.png", width = 7, height = 4, dpi = 1200) }

Gx2x2Fit$coefficients
print(Gx2x2Fit$plot) 
 
if (save_png) { ggsave("Figures/QHO_x2x2_correlator.png", width = 7, height = 4, dpi = 1200) }

# Choose noiseless region based on plateau from log ratio

{
  Gx1x1Fit <- fit_correlator_exp(correlator = Gx1x1Data, noiseless_region = noiseless_region_1, A_guess = 0.5, DeltaE_guess = 1.0, latticeSpacing = latticeSpacing)
  Gx2x2Fit <- fit_correlator_exp(correlator = Gx2x2Data, noiseless_region = noiseless_region_2, A_guess = 0.5, DeltaE_guess = 2.0, latticeSpacing = latticeSpacing)
}

Gx1x1Fit$coefficients
print(Gx1x1Fit$plot)

if (save_png) { ggsave("Figures/DWP_x1x1_correlator.png", width = 7, height = 4, dpi = 1200) }

Gx2x2Fit$coefficients
print(Gx2x2Fit$plot)
 
if (save_png) { ggsave("Figures/DWP_x2x2_correlator.png", width = 7, height = 4, dpi = 1200) }

# nolint end