# nolint start

#####################################################################################################################################
#################################################### Double Well Potential ##########################################################
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
  
  accRateThermData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/accRateTherm"))))
  accRateData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/accRate"))))
  E0ThermData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/E0Therm"))))
  E0Data <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/E0"))))
  histogramData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/histogram"))))
  Gx1x1Data <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/Gx1x1"))))
  Gx2x2Data <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/Gx2x2"))))
  instantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/instantons"))))
  antiInstantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/antiInstantons"))))
  headerData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/headerInfo"))))

  diagEnergiesData <- read.csv("DWP diagonalisation/DWP diagonalisation/energies.csv")
  diagWFData <- read.csv("DWP diagonalisation/DWP diagonalisation/wavefunctions.csv")
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
  wellCentres <- headerData[11]
  lambdaDWP <- headerData[12]

  beta <- pathLength * latticeSpacing
  thermMeasures <- thermSweeps / thermInterval

  # Diagonalisation Energies
  energy_row    <- match(wellCentres, diagEnergiesData[, 1])
  E0_diag       <- diagEnergiesData[energy_row, 2]
  E1_diag       <- diagEnergiesData[energy_row, 3]
  E2_diag       <- diagEnergiesData[energy_row, 4]
  Split_diag    <- diagEnergiesData[energy_row, 5]
  Split_diag_2  <- diagEnergiesData[energy_row, 6]

  # ABCs analysis
  omegaDWP  <- sqrt(8 * lambdaDWP * wellCentres^2)

  S_inst    <- (2/3) * omegaDWP * (wellCentres^2) 
  alpha     <- 1 / 12 # A complicated calculation performed in Zinn-Justin 1993 (or ABCs of Instantons) gives this value
  K         <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5) # A prefactor for the splitting energy
  Split_ABCs <- K * exp(-S_inst)
  E0_ABCs    <- 0.5 * omegaDWP - (Split_ABCs / 2)
  E1_ABCs    <- 0.5 * omegaDWP + (Split_ABCs / 2)
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
    labs(x = "Sweep", y = expression("Ground State Energy" ~ E[0]), title = "Average Thermalisation of the DWP") +
    theme_minimal(base_size = 14)
}

if (save_png) { ggsave("Figures/DWP_thermalisation.png", width = 6, height = 4, dpi = 1200) }

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
        labs(x = "Measure", y = expression("Ground State Energy" ~ E[0]), title = "Average Decorrelation of the DWP") +
        theme_minimal(base_size = 14)
}

if (save_png) { ggsave("Figures/DWP_decorrelation.png", width = 6, height = 4, dpi = 1200) }

###~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Auto-correlation function ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~###

{
  acfResult <- acf(E0Data[1:measures], lag.max = measures - 1, plot = TRUE)
  acfVals <- acfResult$acf
  mean(abs(as.vector(acfVals)[2:length(acfVals)])) # For decorrelated data, we want low acf values (roughly < 0.1)
}

if (save_png) {
    png("Figures/acf_plot_DWP.png", width = 800, height = 600)

    acfResult <- acf(E0Data[1:measures], lag.max = measures - 1, plot = FALSE)

    plot(acfResult, main = "Autocorrelation of Ground State Energy for the DWP", xlab = "Lag (sweeps)", ylab = "Autocorrelation")

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

  E0_pc_error <- abs(((E0 - E0_diag) / E0_diag) * 100)
}

E0; E0_diag; mean(E0RepeatAvg) + E0StandardError; mean(E0RepeatAvg) - E0StandardError ; E0_pc_error

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
      labs(x = expression("Ground State Energy " ~ E[0]), y = "Probability Density") +
      theme_minimal(base_size = 14) +
      theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

    shapiroTest <- shapiro.test(E0RepeatAvg) # Run a Shapiro test on the data, values of p < 0.05 are not normally distributed
    qqPlot <- ggplot(data.frame(E0RepeatAvg), aes(sample = E0RepeatAvg)) +
      stat_qq() +
      stat_qq_line(color = "red") +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_minimal(base_size = 14) +
      theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

    grid.arrange(histPlot, qqPlot, ncol = 2)
}
save_png <- TRUE
if (save_png) {
      png("Figures/DWP_Histogram_qq.png", width = 1200, height = 800)
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
  
    sigmaDWP <- 1 / sqrt(omegaDWP)
    xMax <- ceiling(wellCentres + 4 * sigmaDWP) + 1
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
    psiAnalytical <- exp(-(omegaDWP / 2) * (x_values + wellCentres)^2) + exp(-(omegaDWP / 2) * (x_values - wellCentres)^2) # Approximate by 2 Gaussians
    
    # Normalise the analytical wave function 
    psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) 

    wave_df$psiAnalytical <- psiAnalytical

    f_target <- wellCentres
    diagWFData <- diagWFData %>%
      filter(abs(f - f_target) < 1e-6)
}

# Comparing wave functions

{
  ggplot(wave_df, aes(x = x)) +
    # MCMC wave function
    geom_point(aes(y = psi, color = "MCMC"), size = 0.6) +

    # Diagonalisation wave function
    geom_line(data = diagWFData, aes(x = x, y = psi0, color = "Diagonalised"), linewidth = 0.5, alpha = 0.7) +

    # ABCs approximation
    geom_line(aes(y = psiAnalytical, color = "Double Gaussian"), linewidth = 0.5, linetype = "dashed") +
    
    labs(x = "Position x", y = expression("Wavefunction " ~ psi[0](x))) +

    scale_color_manual(name = "", values = c("MCMC" = "red", "Diagonalised" = "black", "Double Gaussian" = "blue")) + 

    theme_minimal(base_size = 14)
}
save_png <- TRUE   
if (save_png) { ggsave("Figures/DWP_wave_functions.png", width = 12, height = 4, dpi = 1200) }

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
    plateau_width = 100
  )
      
  E1 <- log_ratio_results_1$E
  delta_E <- log_ratio_results_1$DeltaE
  # Estimate error by the standard deviation of ratios in the best plateau 
  delta_E_err <- sd(log_ratio_results_1$plateau) / sqrt(length(log_ratio_results_1$plateau))

  noiseless_region_1 <- log_ratio_results_1$plateau_indices

  E1_pc_error <- abs(((E1 - E1_diag) / E1_diag) * 100)
  delta_E_pc_error <- abs(((delta_E - Split_diag) / Split_diag) * 100)
}

E1; E1_diag; E1 - delta_E_err - E0StandardError; E1 + delta_E_err + E0StandardError; E1_pc_error
delta_E; Split_diag; delta_E - delta_E_err; delta_E + delta_E_err; delta_E_pc_error

# E2

{
  log_ratio_results_2 <- log_ratio_energy(
    correlator = Gx2x2Data,
    latticeSpacing = latticeSpacing,
    E0 = E0,
    plateau_width = 20
  )
      
  E2 <- log_ratio_results_2$E
  delta_E_2 <- log_ratio_results_2$DeltaE
  # Estimate error by the standard deviation of ratios in the best plateau 
  delta_E_2_err <- sd(log_ratio_results_2$plateau) / sqrt(length(log_ratio_results_2$plateau))

  noiseless_region_2 <- log_ratio_results_2$plateau_indices

  E2_pc_error <- abs(((E2 - E2_diag) / E2_diag) * 100)
  delta_E_2_pc_error <- abs(((delta_E_2 - Split_diag_2) / Split_diag_2) * 100)
}

E2; E2_diag; E2 - delta_E_2_err - E0StandardError; E2 + delta_E_2_err + E0StandardError; E2_pc_error
delta_E_2; Split_diag_2; delta_E_2 - delta_E_2_err; delta_E_2 + delta_E_2_err; delta_E_2_pc_error

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
  Gx1x1Fit <- fit_correlator_exp(corr = Gx1x1Data, noiseless_region = 50:500, A_guess = E0_diag, DeltaE_guess = Split_diag, latticeSpacing = latticeSpacing)
  Gx2x2Fit <- fit_correlator_exp(corr = Gx2x2Data, noiseless_region = 1:30, A_guess = E0_diag, DeltaE_guess = Split_diag_2, latticeSpacing = latticeSpacing)
}

Gx1x1Fit$coefficients
print(Gx1x1Fit$plot)

if (save_png) { ggsave("Figures/DWP_x1x1_correlator.png", width = 7, height = 4, dpi = 1200) }

Gx2x2Fit$coefficients
print(Gx2x2Fit$plot)
 
if (save_png) { ggsave("Figures/DWP_x2x2_correlator.png", width = 7, height = 4, dpi = 1200) }

# Choose noiseless region based on plateau from log ratio

{
  Gx1x1Fit <- fit_correlator_exp(corr = Gx1x1Data, noiseless_region = noiseless_region_1, A_guess = E0_diag, DeltaE_guess = Split_diag, latticeSpacing = latticeSpacing)
  Gx2x2Fit <- fit_correlator_exp(corr = Gx2x2Data, noiseless_region = noiseless_region_2, A_guess = E0_diag, DeltaE_guess = Split_diag_2, latticeSpacing = latticeSpacing)
}

Gx1x1Fit$coefficients
print(Gx1x1Fit$plot)

if (save_png) { ggsave("Figures/DWP_x1x1_correlator.png", width = 7, height = 4, dpi = 1200) }

Gx2x2Fit$coefficients
print(Gx2x2Fit$plot)
 
if (save_png) { ggsave("Figures/DWP_x2x2_correlator.png", width = 7, height = 4, dpi = 1200) }

###~~~~~~~~~~~~~~~~###
### Tunneling data ###
###~~~~~~~~~~~~~~~~###

# Instanton density

mean(instantonsData) / pathLength

# An approximation for the instanton action based on the number of tunnelling events  

{
  S_inst_approximation <- log(pathLength / mean(instantonsData)) 
  S_inst_approximation <- S_inst_approximation * omegaDWP / pi * sqrt(S_inst_approximation)
} 

S_inst; S_inst_approximation

# Comparing energies against diagonalisation method and ABCs

E0; E0_diag; E0_ABCs

E1; E1_diag; E1_ABCs

E1 - E0; Split_diag; Split_ABCs

#####################################################################################################################################
################################################# Vary Well Centres and Beta ########################################################
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

  read_h5 <- function(name) as.numeric(unlist(h5read(dataFile, paste0(base, name))))
  
  well_centres_vec <- c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
  betas_vec <- c(500, 400, 250, 100, 75, 50, 25)
  latticeSpacing <- 0.05 
  path_lengths <- as.integer(betas_vec / latticeSpacing)
  
  wellReps <- length(well_centres_vec)
  betaReps <- length(betas_vec)
  
  # Build paths matrix (paths are always same size, so matrix is fine here)
  full_paths <- matrix("", nrow = wellReps, ncol = betaReps)
  
  for (i in seq_along(well_centres_vec)) {
    wc_str <- sprintf("%.6f", well_centres_vec[i])
    for (j in seq_along(betas_vec)) {
      ls_str <- sprintf("%.6f", betas_vec[j])
      full_paths[i, j] <- paste0("/DWP/six/", wc_str, "/", ls_str, "/")
    }
  }
  
  # Allocate data as nested lists
  accRateData        <- vector("list", wellReps)
  E0Data             <- vector("list", wellReps)
  histogramData      <- vector("list", wellReps)
  Gx1x1Data          <- vector("list", wellReps)
  Gx2x2Data          <- vector("list", wellReps)
  instantonsData     <- vector("list", wellReps)
  antiInstantonsData <- vector("list", wellReps)
  headerData         <- vector("list", wellReps)
  
  # Read data
  for (i in seq_along(well_centres_vec)) {
    
    # Initialize sublists for each well center
    accRateData[[i]]        <- vector("list", betaReps)
    E0Data[[i]]             <- vector("list", betaReps)
    histogramData[[i]]      <- vector("list", betaReps)
    Gx1x1Data[[i]]          <- vector("list", betaReps)
    Gx2x2Data[[i]]          <- vector("list", betaReps)
    instantonsData[[i]]     <- vector("list", betaReps)
    antiInstantonsData[[i]] <- vector("list", betaReps)
    headerData[[i]]         <- vector("list", betaReps)
    
    for (j in seq_along(betas_vec)) {
      base <- full_paths[i, j]
      
      # Read everything into the nested list structure
      accRateData[[i]][[j]]        <- read_h5("accRate")
      E0Data[[i]][[j]]             <- read_h5("E0")
      histogramData[[i]][[j]]      <- read_h5("histogram")
      Gx1x1Data[[i]][[j]]          <- read_h5("Gx1x1")
      Gx2x2Data[[i]][[j]]          <- read_h5("Gx2x2")
      instantonsData[[i]][[j]]     <- read_h5("instantons")
      antiInstantonsData[[i]][[j]] <- read_h5("antiInstantons")
      headerData[[i]][[j]]         <- read_h5("headerInfo")
    }
  }
  # Data is now stored as Variable[[well_index]][[beta_index]]
}

###~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Theoretical calculations ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~###

{
  accRateInterval <- headerData[[1]][[1]][4]
  decorrSweeps <- headerData[[1]][[1]][5]
  thermSweeps <- headerData[[1]][[1]][6]
  thermInterval <- headerData[[1]][[1]][7]
  measures <- headerData[[1]][[1]][8]
  repeats <- headerData[[1]][[1]][9]
  numBins <- headerData[[1]][[1]][10]
  lambdaDWP <- headerData[[1]][[1]][12]

  wellCentres <- well_centres_vec

  omegaDWP <- sqrt(8 * lambdaDWP * wellCentres^2)
  S_inst   <- (2/3) * omegaDWP * (wellCentres^2)

  alpha <- 1 / 12
  K <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5)

  Split_ABCs <- K * exp(-S_inst)
  E0_ABCs <- 0.5 * omegaDWP - (Split_ABCs / 2)
  E1_ABCs <- 0.5 * omegaDWP + (Split_ABCs / 2)

  # Grabovsky Splitting energy
  Split_Grabovsky <- (omegaDWP / pi) * exp(-omegaDWP * wellCentres^2)
  E0_Grabovsky    <- 0.5 * omegaDWP - (Split_Grabovsky / 2)
  E1_Grabovsky    <- 0.5 * omegaDWP + (Split_Grabovsky / 2)

  # Diagonalisation data
  diagEnergiesData <- read.csv("DWP diagonalisation/DWP diagonalisation/energies.csv")

  energy_row <- match(wellCentres, diagEnergiesData[, 1])
  E0_diag       <- diagEnergiesData[energy_row, 2]
  E1_diag       <- diagEnergiesData[energy_row, 3]
  E2_diag       <- diagEnergiesData[energy_row, 4]
  Split_diag    <- diagEnergiesData[energy_row, 5]
  Split_diag_2 <- diagEnergiesData[energy_row, 6]
}

###~~~~~~~~~~~~~~~~~###
### Acceptance rate ###
###~~~~~~~~~~~~~~~~~###

{
  acc_matrix <- matrix(0, nrow = wellReps, ncol = betaReps)

  for (i in 1:wellReps) {
    for (j in 1:betaReps) {
      acc_matrix[i, j] <- round(mean(accRateData[[i]][[j]]) * 100, digits = 1)
    }
  }

  rownames(acc_matrix) <- paste0(well_centres_vec)
  colnames(acc_matrix) <- paste0(betas_vec)
  print(acc_matrix)
}

###~~~~~~~~~~~~~~~~~~~~~~~~###
### Calculate ground state ### 
###~~~~~~~~~~~~~~~~~~~~~~~~###

{
  E0_vals   <- matrix(0, nrow = wellReps, ncol = betaReps)
  E0_errors <- matrix(0, nrow = wellReps, ncol = betaReps)

  for (i in seq_len(wellReps)) {
    for (j in seq_len(betaReps)) {

      E0_series <- E0Data[[i]][[j]]

      E0Split <- split(E0_series, rep(1:repeats, each = measures))
      
      E0RepeatAvg <- sapply(E0Split, mean)
      
      E0_vals[i, j]   <- mean(E0RepeatAvg)
      E0_errors[i, j] <- sd(E0RepeatAvg) / sqrt(length(E0RepeatAvg))
    }
  }

  rownames(E0_vals)   <- well_centres_vec
  colnames(E0_vals)   <- betas_vec
  rownames(E0_errors) <- well_centres_vec
  colnames(E0_errors) <- betas_vec
}

# Ground state energies
print(E0_vals, digits = 4)

{
  plot_df <- expand.grid(
    WellCentre = well_centres_vec,
    Beta = betas_vec
  )

  plot_df$Energy <- as.vector(E0_vals)
  plot_df$Error  <- as.vector(E0_errors)

  plot_df <- plot_df %>%
    mutate(
      # Match the diagonalised ground state energy for each well centre
      E0_diag = E0_diag[match(WellCentre, well_centres_vec)],
      Well_Color = as.factor(WellCentre)
    )
}

error_summary <- plot_df %>%
  filter(Beta == 250) %>% 
  mutate(
    # Statistical Percentage Error (MCMC noise)
    stat_pc_error = (Error / Energy) * 100,
    
    # Systematic Percentage Bias (Diagonalization used as benchmark)
    sys_pc_error = (abs(Energy - E0_diag) / E0_diag) * 100
  )

# Calculate the averages across all well centres
avg_stat_error <- mean(error_summary$stat_pc_error)
avg_sys_error  <- mean(error_summary$sys_pc_error)

# Print results
cat(sprintf("Average Statistical Error: %.3f%%\n", avg_stat_error))
cat(sprintf("Average Systematic Bias:    %.3f%%\n", avg_sys_error))

# Combine in quadrature
avg_total_error <- sqrt(avg_stat_error^2 + avg_sys_error^2)

cat(sprintf("Total Combined Relative Error: %.3f%%\n", avg_total_error))

ggplot(plot_df, aes(x = Beta, y = Energy, color = Well_Color)) +
  # Diagonalised values
  geom_hline(aes(yintercept = E0_diag, color = Well_Color), linetype = "dashed") +
  
  # MCMC ground state energy
  geom_line(aes(group = WellCentre), linewidth = 0.7) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Energy - Error, ymax = Energy + Error), width = 0.03) +
  
  # Apply log scale to x-axis as Beta spans orders of magnitude
  scale_x_log10(breaks = betas_vec) +
  
  # Labels and theme
  labs(title = "Ground State Convergence with Decreasing Temperature", x = expression("Inverse Temperature"~Beta), 
      y = expression("Ground State Energy " ~ E[0]), color = "Well Centre") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

if (save_png) { ggsave("Figures/DWP_E0_vary_beta.png", width = 11, height = 10, dpi = 1200) }

###~~~~~~~~~~~~~~~~~~~~~~~###
### Excited energy states ###
###~~~~~~~~~~~~~~~~~~~~~~~###

# Log ratio method

find_plateau <- function(ratios, width = 10) {
  n <- length(ratios)
  
  # Check if data is shorter than width
  if (n <= width) return(1:n)
  
  # Calculate variance for every possible window of size 'width'
  # (n - width + 1) is the number of possible windows
  variances <- sapply(1:(n - width + 1), function(i) {
    var(ratios[i:(i + width - 1)])
  })
  
  # Find the start index of the window with the minimum variance
  best_start <- which.min(variances)
  
  # Return the full range of indices for that plateau
  return(best_start:(best_start + width - 1))
}

log_ratio_energy <- function(correlator, latticeSpacing, E0, plateau_width = 10) {
  # Calculate all available log-ratios
  max_lag <- length(correlator) / 2
  raw_ratios <- numeric(max_lag)
  for (t in 1:max_lag) {
    if (correlator[t] > 0 && correlator[t + 1] > 0) {
      raw_ratios[t] <- log(correlator[t] / correlator[t + 1]) / latticeSpacing
    } else {
      raw_ratios[t] <- NA
    }
  }
  
  clean_ratios <- raw_ratios[!is.na(raw_ratios)]
  
  # Find the flattest window of the specified width
  plateau_indices <- find_plateau(clean_ratios, width = plateau_width)
  plateau_values  <- clean_ratios[plateau_indices]
  
  # Compute results
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

{
  E1_vals_log       <- matrix(NA, wellReps, betaReps)
  delta_E_vals_log  <- matrix(NA, wellReps, betaReps)
  delta_E_err_log  <- matrix(NA, wellReps, betaReps)

  noiseless_regions <- matrix(NA, wellReps, betaReps)

  for (i in seq_len(wellReps)) {
    for (j in seq_len(betaReps)) {
      # Log-ratio estimation for E1
      log_ratio_results <- log_ratio_energy(
        correlator = Gx1x1Data[[i]][[j]],
        latticeSpacing = latticeSpacing,
        E0 = E0_vals[i, j],
        plateau_width = 5
      )
      
      E1_vals_log[i, j] <- log_ratio_results$E
      delta_E_vals_log[i, j] <- log_ratio_results$DeltaE
      # Estimate error by the standard deviation of ratios in the best plateau 
      delta_E_err_log[i, j]  <- sd(log_ratio_results$plateau) / sqrt(length(log_ratio_results$plateau))

      # noiseless_regions[i, j] <- log_ratio_results$plateau_indices
    }
  }
}

delta_E_vals_log

{
  plot_log_df <- expand.grid(
    WellCentre = well_centres_vec,
    Beta = betas_vec
  )

  plot_log_df$E1              <- as.vector(E1_vals_log)
  plot_log_df$Splitting       <- as.vector(delta_E_vals_log)
  plot_log_df$Splitting_err   <- as.vector(delta_E_err_log)

  plot_log_df <- plot_log_df %>%
    mutate(
      Well_Factor = as.factor(WellCentre)
    )
}

ggplot(plot_log_df, aes(x = Beta, y = Splitting, color = Well_Factor)) +
  # MCMC data
  geom_line(aes(group = Well_Factor), size = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = Splitting - Splitting_err, ymax = Splitting + Splitting_err), width = 0.01) +
  
  # Diagonalisation data
  geom_hline(aes(yintercept = Split_diag[Well_Factor], color = Well_Factor), linetype = "dashed") +
  
  # Apply log scale to x-axis as Beta spans orders of magnitude
  scale_x_log10(breaks = betas_vec) +

  labs(
    title = "Energy Splitting Convergence",
    x = expression("Log of Inverse Temperature" ~ Log[10](Beta)),
    y = expression(Delta * E),
    color = "Well Centre") +
  # ylim(0, 1) + 
    theme_minimal(base_size = 14)

if (save_png) { ggsave("Figures/DWP_delta_E_vary_beta.png", width = 10, height = 10, dpi = 1200) }

# Exponential fit method

fit_correlator <- function(correlator, latticeSpacing, noiseless_region, E0, DeltaE_guess) {
    
  correlator <- correlator[noiseless_region] # Use only noiseless region of the correlator
  
  df <- data.frame(lag = noiseless_region, correlation = correlator)
  
  fit <- nlsLM(
    correlation ~ A * exp(-DeltaE * lag * latticeSpacing) + c,
    data = df,
    start = list(A = correlator[1], c = min(correlator), DeltaE = DeltaE_guess),
    control = nls.lm.control(maxiter = 1024)
  )
  
  list(
    E = E0 + coef(fit)["DeltaE"],
    DeltaE = coef(fit)["DeltaE"],
    fit = predict(fit, newdata = df),
    df = df
  )
}

plot_correlator <- function(correlator, noiseless_region, latticeSpacing, wellCentres, exponential_fit_results = NULL) {
  n_total <- length(correlator)
  beta <- latticeSpacing * n_total
  
  full_time_axis <- (0:(n_total - 1)) * latticeSpacing
  
  df <- data.frame(
    time = full_time_axis[noiseless_region], 
    correlator = correlator[noiseless_region]
  )
  
  p <- ggplot(df, aes(x = time, y = correlator)) +
    geom_line(linewidth = 0.5) + 
    geom_point(size = 1) +
    labs(x = "Euclidean time t", y = "G(t)", 
          title = bquote("DWP correlator for " ~ Beta == .(beta) ~ 
                        ", a =" ~ .(latticeSpacing) ~ 
                        ", wellCentres =" ~ .(wellCentres))) +
    theme_minimal(base_size = 14)
  
  if (!is.null(exponential_fit_results)) {
    df_fit <- data.frame(
      time = full_time_axis[noiseless_region],
      fit_val = exponential_fit_results$fit
    )
    
    p <- p + geom_line(data = df_fit, aes(x = time, y = fit_val), 
                        color = "red", linetype = "dashed")
  }
  
  return(p)
}

{
  E1_vals_fit       <- matrix(NA, wellReps, betaReps)
  delta_E_vals_fit  <- matrix(NA, wellReps, betaReps)
  delta_E_err_fit   <- matrix(NA, wellReps, betaReps)
}

# Naively choose noiseless region

{
  naive_noiseless_region <- 1:20 

  for (i in seq_len(wellReps)) {
    for (j in seq_len(betaReps)) {

      # Nonlinear Fit for Energy Splitting (Delta E)
      exponential_fit_results <- fit_correlator(
        Gx1x1Data[[i]][[j]],
        latticeSpacing,
        naive_noiseless_region,
        E0_vals[i, j],
        DeltaE_guess = Split_diag[i] # Could estimate splitting energy if no known values were available
      )
      
      E1_vals_fit[i, j] <- exponential_fit_results$E
      delta_E_vals_fit[i, j] <- exponential_fit_results$DeltaE
      delta_E_err_fit[i, j]  <- sd(exponential_fit_results$DeltaE) / sqrt(length(exponential_fit_results$DeltaE))
    }
  }
}

delta_E_vals_fit

# Choose noiseless region based on plateau from log ratio

{
  for (i in seq_len(wellReps)) {
    for (j in seq_len(betaReps)) {

      # Nonlinear Fit for Energy Splitting (Delta E)
      exponential_fit_results <- fit_correlator(
        Gx1x1Data[[i]][[j]],
        latticeSpacing,
        noiseless_regions[i, j],
        E0_vals[i, j],
        DeltaE_guess = Split_diag[i] # Could estimate splitting energy if no known values were available
      )
      
      E1_vals_fit[i, j] <- exponential_fit_results$E
      delta_E_vals_fit[i, j] <- exponential_fit_results$DeltaE
      delta_E_err_fit[i, j]  <- sd(exponential_fit_results$DeltaE) / sqrt(length(exponential_fit_results$DeltaE))
    }
  }
}

delta_E_vals_fit

# Plot any of the correlator fits to check the fit
wellIndex <- 1
betaIndex <- 1

plot_correlator(Gx1x1Data[[wellIndex]][[betaIndex]], naive_noiseless_region, latticeSpacing, well_centres_vec[wellIndex],
  fit_correlator(
      Gx1x1Data[[wellIndex]][[betaIndex]],
      latticeSpacing,
      naive_noiseless_region,
      E0_vals[wellIndex, betaIndex],
      DeltaE_guess = Split_diag[wellIndex]))

plot_correlator(Gx1x1Data[[wellIndex]][[betaIndex]], noiseless_regions[wellIndex, betaIndex], latticeSpacing, well_centres_vec[wellIndex],
  fit_correlator(
      Gx1x1Data[[wellIndex]][[betaIndex]],
      latticeSpacing,
      noiseless_regions[wellIndex, betaIndex],
      E0_vals[wellIndex, betaIndex],
      DeltaE_guess = Split_diag[wellIndex]))

#####################################################################################################################################
########################################## Vary Well Centres and Lattice Spacing ####################################################
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

###~~~~~~~~~~~###
### Read data ###
###~~~~~~~~~~~###

{
  dataFile <- "data2.h5"
  
  read_h5 <- function(name) as.numeric(unlist(h5read(dataFile, paste0(base, name))))

  well_centres_vec <- c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
  # Now we vary lattice spacing (a) instead of beta
  lattice_spacing_vec <- c(0.5, 0.4, 0.3, 0.25, 0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05)
  beta_val <- 500 # Keep beta constant for this comparison
  
  wellReps <- length(well_centres_vec)
  lsReps   <- length(lattice_spacing_vec)
  
  # Build paths matrix indexed by [WellCentre, LatticeSpacing]
  full_paths <- matrix("", nrow = wellReps, ncol = lsReps)
  
  for (i in seq_along(well_centres_vec)) {
    wc_str <- sprintf("%.6f", well_centres_vec[i])
    for (j in seq_along(lattice_spacing_vec)) {
      ls_str <- sprintf("%.6f", lattice_spacing_vec[j])
      full_paths[i, j] <- paste0("/DWP/seven/", wc_str, "/", ls_str, "/")
    }
  }
  
  # Allocate lists for data
  accRateData        <- vector("list", wellReps)
  E0Data             <- vector("list", wellReps)
  histogramData      <- vector("list", wellReps)
  Gx1x1Data          <- vector("list", wellReps)
  Gx2x2Data          <- vector("list", wellReps)
  instantonsData     <- vector("list", wellReps)
  antiInstantonsData <- vector("list", wellReps)
  headerData         <- vector("list", wellReps)
  
  # Read data
  for (i in seq_along(well_centres_vec)) {
    
    # Initialize sublists for each well center
    accRateData[[i]]        <- vector("list", lsReps)
    E0Data[[i]]             <- vector("list", lsReps)
    histogramData[[i]]      <- vector("list", lsReps)
    Gx1x1Data[[i]]          <- vector("list", lsReps)
    Gx2x2Data[[i]]          <- vector("list", lsReps)
    instantonsData[[i]]     <- vector("list", lsReps)
    antiInstantonsData[[i]] <- vector("list", lsReps)
    headerData[[i]]         <- vector("list", lsReps)
    
    for (j in seq_along(lattice_spacing_vec)) {
      base <- full_paths[i, j]
      
      # Read everything into the nested list structure
      accRateData[[i]][[j]]        <- read_h5("accRate")
      E0Data[[i]][[j]]             <- read_h5("E0")
      histogramData[[i]][[j]]      <- read_h5("histogram")
      Gx1x1Data[[i]][[j]]          <- read_h5("Gx1x1")
      Gx2x2Data[[i]][[j]]          <- read_h5("Gx2x2")
      instantonsData[[i]][[j]]     <- read_h5("instantons")
      antiInstantonsData[[i]][[j]] <- read_h5("antiInstantons")
      headerData[[i]][[j]]         <- read_h5("headerInfo")
    }
  }
}

###~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Theoretical calculations ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~###

{
  accRateInterval <- headerData[[1]][[1]][4]
  decorrSweeps <- headerData[[1]][[1]][5]
  thermSweeps <- headerData[[1]][[1]][6]
  thermInterval <- headerData[[1]][[1]][7]
  measures <- headerData[[1]][[1]][8]
  repeats <- headerData[[1]][[1]][9]
  numBins <- headerData[[1]][[1]][10]
  lambdaDWP <- headerData[[1]][[1]][12]

  wellCentres <- well_centres_vec

  omegaDWP <- sqrt(8 * lambdaDWP * wellCentres^2)
  S_inst   <- (2/3) * omegaDWP * (wellCentres^2)

  alpha <- 1 / 12
  K <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5)

  Split_ABCs <- K * exp(-S_inst)
  E0_ABCs <- 0.5 * omegaDWP - (Split_ABCs / 2)
  E1_ABCs <- 0.5 * omegaDWP + (Split_ABCs / 2)

  # Diagonalisation data
  diagEnergiesData <- read.csv("DWP diagonalisation/DWP diagonalisation/energies.csv")

  energy_row <- match(wellCentres, diagEnergiesData[, 1])
  E0_diag       <- diagEnergiesData[energy_row, 2]
  E1_diag       <- diagEnergiesData[energy_row, 3]
  E2_diag       <- diagEnergiesData[energy_row, 4]
  Split_diag    <- diagEnergiesData[energy_row, 5]
  Split_diag_2 <- diagEnergiesData[energy_row, 6]
}

###~~~~~~~~~~~~~~~~~~~~~~~~###
### Calculate ground state ### 
###~~~~~~~~~~~~~~~~~~~~~~~~###

{
  E0_vals   <- matrix(0, nrow = wellReps, ncol = lsReps)
  E0_errors <- matrix(0, nrow = wellReps, ncol = lsReps)

  for (i in seq_len(wellReps)) {
    for (j in seq_len(lsReps)) {

      E0_series <- E0Data[[i]][[j]]

      E0Split <- split(E0_series, rep(1:repeats, each = measures))
      
      E0RepeatAvg <- sapply(E0Split, mean)
      
      E0_vals[i, j]   <- mean(E0RepeatAvg)
      E0_errors[i, j] <- sd(E0RepeatAvg) / sqrt(length(E0RepeatAvg))
    }
  }

  rownames(E0_vals)   <- well_centres_vec
  colnames(E0_vals)   <- lattice_spacing_vec
  rownames(E0_errors) <- well_centres_vec
  colnames(E0_errors) <- lattice_spacing_vec
}

# Ground state energies
print(E0_vals, digits = 3)

{
  plot_df <- expand.grid(
    WellCentre = well_centres_vec,
    latticeSpacing = lattice_spacing_vec
  )

  plot_df$Energy <- as.vector(E0_vals)
  plot_df$Error  <- as.vector(E0_errors)

  plot_df <- plot_df %>%
    mutate(
      # Match the diagonalised ground state energy for each well centre
      E0_diag = E0_diag[match(WellCentre, well_centres_vec)],
      Well_Color = as.factor(WellCentre)
    )
}

ls_error_summary <- plot_df %>%
  filter(latticeSpacing == 0.075) %>% 
  mutate(
    stat_pc_error = (Error / Energy) * 100,
    sys_pc_error  = (abs(Energy - E0_diag) / E0_diag) * 100,
    total_pc_err  = sqrt(stat_pc_error^2 + sys_pc_error^2)
  )

# If you want the average across ALL lattice spacings tested:
avg_ls_stat <- mean(ls_error_summary$stat_pc_error)
avg_ls_sys  <- mean(ls_error_summary$sys_pc_error)
avg_ls_tot  <- mean(ls_error_summary$total_pc_err)

# Print results
cat(sprintf("Average Stats Error (all a): %.3f%%\n", avg_ls_stat))
cat(sprintf("Average Sys Bias (all a):    %.3f%%\n", avg_ls_sys))
cat(sprintf("Total Combined Error:        %.3f%%\n", avg_ls_tot))

ggplot(plot_df, aes(x = 1 / (latticeSpacing^2), y = Energy, color = Well_Color)) +
  # Diagonalised values
  geom_hline(aes(yintercept = E0_diag, color = Well_Color), linetype = "dashed") +
  
  # MCMC ground state energy
  geom_line(aes(group = WellCentre), size = 0.7) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Energy - Error, ymax = Energy + Error), width = 10) +
  
  # Labels and theme
  labs(title = "Ground State Convergence with Decreasing Lattice Spacing", x = expression("Inverse square lattice spacing" ~ 1 / a^2), 
      y = expression("Ground State Energy " ~ E[0]), color = "Well Centre") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

if (save_png) { ggsave("Figures/DWP_E0_vary_ls.png", width = 11, height = 10, dpi = 1200) }

###~~~~~~~~~~~~~~~~~~~~~~~###
### Excited energy states ###
###~~~~~~~~~~~~~~~~~~~~~~~###

# Log-Ratio Method

find_plateau <- function(ratios, width = 10) {
  n <- length(ratios)
  
  # Check if data is shorter than width
  if (n <= width) return(1:n)
  
  # Calculate variance for every possible window of size 'width'
  # (n - width + 1) is the number of possible windows
  variances <- sapply(1:(n - width + 1), function(i) {
    var(ratios[i:(i + width - 1)])
  })
  
  # Find the start index of the window with the minimum variance
  best_start <- which.min(variances)
  
  # Return the full range of indices for that plateau
  return(best_start:(best_start + width - 1))
}

log_ratio_energy <- function(correlator, latticeSpacing, E0, plateau_width = 10) {
  # Calculate all available log-ratios
  max_lag <- length(correlator) / 2
  raw_ratios <- numeric(max_lag)
  for (t in 1:max_lag) {
    if (correlator[t] > 0 && correlator[t + 1] > 0) {
      raw_ratios[t] <- log(correlator[t] / correlator[t + 1]) / latticeSpacing
    } else {
      raw_ratios[t] <- NA
    }
  }
  
  clean_ratios <- raw_ratios[!is.na(raw_ratios)]
  
  # Find the flattest window of the specified width
  plateau_indices <- find_plateau(clean_ratios, width = plateau_width)
  plateau_values  <- clean_ratios[plateau_indices]
  
  # Compute results
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

{
  E1_vals_log       <- matrix(NA, wellReps, lsReps)
  delta_E_vals_log  <- matrix(NA, wellReps, lsReps)
  delta_E_err_log  <- matrix(NA, wellReps, lsReps)

  noiseless_regions <- matrix(NA, wellReps, lsReps)

  for (i in seq_len(wellReps)) {
    for (j in seq_len(lsReps)) {
      # Log-ratio estimation for E1
      log_ratio_results <- log_ratio_energy(
        correlator = Gx1x1Data[[i]][[j]],
        latticeSpacing = lattice_spacing_vec[[j]],
        E0 = E0_vals[i, j],
        plateau_width = 10
      )
      
      E1_vals_log[i, j] <- log_ratio_results$E
      delta_E_vals_log[i, j] <- log_ratio_results$DeltaE
      # Estimate error by the standard deviation of ratios in the best plateau 
      delta_E_err_log[i, j]  <- sd(log_ratio_results$plateau) / sqrt(length(log_ratio_results$plateau))

      # noiseless_regions[i, j] <- log_ratio_results$plateau_indices
    }
  }
}

delta_E_vals_log

{
  plot_log_df <- expand.grid(
    WellCentre = well_centres_vec,
    LatticeSpacing = lattice_spacing_vec
  )

  plot_log_df$E1              <- as.vector(E1_vals_log)
  plot_log_df$Splitting       <- as.vector(delta_E_vals_log)
  plot_log_df$Splitting_err   <- as.vector(delta_E_err_log)

  delta_E_err_log

  plot_log_df <- plot_log_df %>%
    mutate(
      Well_Factor = as.factor(WellCentre)
    )
}

ggplot(plot_log_df, aes(x = 1 / (LatticeSpacing^2), y = Splitting, color = Well_Factor)) +
  # MCMC data
  geom_line(aes(group = Well_Factor), size = 0.8) +
  geom_point(size = 2.5) +
  # geom_errorbar(aes(ymin = Splitting - Splitting_err, ymax = Splitting + Splitting_err), width = 0.5) +
  
  # Diagonalisation data
  geom_hline(aes(yintercept = Split_diag[Well_Factor], color = Well_Factor), linetype = "dashed") +

  labs(
    title = "Energy Splitting Convergence against Lattice Spacing",
    x = expression("Inverse square lattice spacing" ~ 1 / (a^2)),
    y = expression(Delta * E),
    color = "Well Centre") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

if (save_png) { ggsave("Figures/DWP_delta_E_vary_ls.png", width = 11, height = 10, dpi = 1200) }

# Exponential fit

fit_correlator <- function(correlator, latticeSpacing, noiseless_region, E0, DeltaE_guess) {
    
  correlator <- correlator[noiseless_region] # Use only noiseless region of the correlator
  
  df <- data.frame(lag = noiseless_region, correlation = correlator)
  
  fit <- nlsLM(
    correlation ~ A * exp(-DeltaE * lag * latticeSpacing) + c,
    data = df,
    start = list(A = correlator[1], c = min(correlator), DeltaE = DeltaE_guess),
    control = nls.lm.control(maxiter = 1024)
  )
  
  list(
    E = E0 + coef(fit)["DeltaE"],
    DeltaE = coef(fit)["DeltaE"],
    fit = predict(fit, newdata = df),
    df = df
  )
}

plot_correlator <- function(correlator, noiseless_region, latticeSpacing, wellCentres, exponential_fit_results = NULL) {
  n_total <- length(correlator)
  beta <- latticeSpacing * n_total
  
  full_time_axis <- (0:(n_total - 1)) * latticeSpacing
  
  df <- data.frame(
    time = full_time_axis[noiseless_region], 
    correlator = correlator[noiseless_region]
  )
  
  p <- ggplot(df, aes(x = time, y = correlator)) +
    geom_line(linewidth = 0.5) + 
    geom_point(size = 1) +
    labs(x = "Euclidean time t", y = "G(t)", 
          title = bquote("DWP correlator for " ~ Beta == .(beta) ~ 
                        ", a =" ~ .(latticeSpacing) ~ 
                        ", wellCentres =" ~ .(wellCentres))) +
    theme_minimal(base_size = 14)
  
  if (!is.null(exponential_fit_results)) {
    df_fit <- data.frame(
      time = full_time_axis[noiseless_region],
      fit_val = exponential_fit_results$fit
    )
    
    p <- p + geom_line(data = df_fit, aes(x = time, y = fit_val), 
                        color = "red", linetype = "dashed")
  }
  
  return(p)
}

{
  E1_vals_fit       <- matrix(NA, wellReps, lsReps)
  delta_E_err_fit   <- matrix(NA, wellReps, lsReps)
  delta_E_vals_fit  <- matrix(NA, wellReps, lsReps)
}

# Naively choose noiseless region

{
  noiseless_region <- 1:150

  for (i in seq_len(wellReps)) {
    for (j in seq_len(lsReps)) {
      if(is.null(Gx1x1Data[[i]][[j]])) next
      
      # Use the specific lattice spacing for this iteration
      current_a <- lattice_spacing_vec[j]
      
      fit_res <- fit_correlator(
        Gx1x1Data[[i]][[j]],
        current_a,
        noiseless_region,
        E0_vals[i, j],
        DeltaE_guess = Split_diag[i] 
      )
      E1_vals_fit[i, j] <- fit_res$E
      delta_E_vals_fit[i, j] <- fit_res$DeltaE
      delta_E_err_fit[i, j] <- sd(fit_res$DeltaE) / sqrt(length(fit_res$DeltaE))
    }
  }
}

delta_E_vals_fit

# Choose noiseless region based on plateau from log ratio

{
  for (i in seq_len(wellReps)) {
    for (j in seq_len(lsReps)) {

      # Nonlinear Fit for Energy Splitting (Delta E)
      exponential_fit_results <- fit_correlator(
        Gx1x1Data[[i]][[j]],
        lattice_spacing_vec[[j]],
        noiseless_regions[i, j],
        E0_vals[i, j],
        DeltaE_guess = Split_diag[i] # Could estimate splitting energy if no known values were available
      )
      
      E1_vals_fit[i, j] <- exponential_fit_results$E
      delta_E_vals_fit[i, j] <- exponential_fit_results$DeltaE
      delta_E_err_fit[i, j]  <- sd(exponential_fit_results$DeltaE) / sqrt(length(exponential_fit_results$DeltaE))
    }
  }
}

delta_E_vals_fit

# Plot any of the correlator fits to check the fit
wellIndex <- 1
lsIndex <- 1

plot_correlator(Gx1x1Data[[wellIndex]][[lsIndex]], noiseless_regions[wellIndex, lsIndex], latticeSpacing, well_centres_vec[wellIndex],
  fit_correlator(
      Gx1x1Data[[wellIndex]][[lsIndex]],
      latticeSpacing,
      noiseless_regions[wellIndex, lsIndex],
      E0_vals[wellIndex, lsIndex],
      DeltaE_guess = Split_diag[wellIndex]))

#####################################################################################################################################
################################################### Vary Well Centres ###############################################################
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

###~~~~~~~~~~~###
### Read data ###
###~~~~~~~~~~~###

{
  dataFile <- "data.h5"
  
  read_h5 <- function(base_path, name) { as.numeric(unlist(h5read(dataFile, paste0(base_path, name)))) }

  well_centres_vec <- c(0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 
                        0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500, 
                        0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725, 0.750, 
                        0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 0.950, 0.975, 1.000, 
                        1.025, 1.050, 1.075, 1.100, 1.125, 1.150, 1.175, 1.200, 1.225, 1.250,  
                        1.275, 1.300, 1.325, 1.350, 1.375, 1.400, 1.425, 1.450, 1.475, 1.500)
  
  wellReps <- length(well_centres_vec)
  
  wc_strs <- sprintf("%.6f", well_centres_vec)
  full_paths <- paste0("/DWP/eight/", wc_strs, "/")
  
  # Allocate lists
  accRateData        <- vector("list", wellReps)
  E0Data             <- vector("list", wellReps)
  histogramData      <- vector("list", wellReps)
  Gx1x1Data          <- vector("list", wellReps)
  Gx2x2Data          <- vector("list", wellReps)
  instantonsData     <- vector("list", wellReps)
  antiInstantonsData <- vector("list", wellReps)
  headerData         <- vector("list", wellReps)
  
  # Read data
  for (i in seq_along(full_paths)) {
    # Extract the specific base path for this iteration
    current_base <- full_paths[i]
    
    accRateData[[i]]        <- read_h5(current_base, "accRate")
    E0Data[[i]]             <- read_h5(current_base, "E0")
    histogramData[[i]]      <- read_h5(current_base, "histogram")
    Gx1x1Data[[i]]          <- read_h5(current_base, "Gx1x1")
    Gx2x2Data[[i]]          <- read_h5(current_base, "Gx2x2")
    instantonsData[[i]]     <- read_h5(current_base, "instantons")
    antiInstantonsData[[i]] <- read_h5(current_base, "antiInstantons")
    headerData[[i]]         <- read_h5(current_base, "headerInfo")
  }
}

###~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Theoretical calculations ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~###

{
  pathLength <- headerData[[1]][1]
  latticeSpacing <- headerData[[1]][2]
  epsilon <- headerData[[1]][3]
  accRateInterval <- headerData[[1]][4]
  # Thermalisation and decorrelation vary but we don't need them
  thermInterval <- headerData[[1]][7]
  measures <- headerData[[1]][8]
  repeats <- headerData[[1]][9]
  numBins <- headerData[[1]][10]
  lambdaDWP <- headerData[[1]][12]

  wellCentres <- well_centres_vec

  # ABCs analysis
  omegaDWP <- sqrt(8 * lambdaDWP * wellCentres^2)
  S_inst   <- (2/3) * omegaDWP * (wellCentres^2)

  alpha <- 1 / 12
  K <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5)

  Split_ABCs <- K * exp(-S_inst)
  E0_ABCs <- 0.5 * omegaDWP - (Split_ABCs / 2)
  E1_ABCs <- 0.5 * omegaDWP + (Split_ABCs / 2)

  # Grabovsky
  Split_Grabovsky <- (omegaDWP / pi) * exp(-omegaDWP * wellCentres^2)
  E0_Grabovsky    <- 0.5 * omegaDWP - (Split_Grabovsky / 2)
  E1_Grabovsky    <- 0.5 * omegaDWP + (Split_Grabovsky / 2)

  # Diagonalisation data
  diagEnergiesData <- read.csv("DWP diagonalisation/DWP diagonalisation/energies.csv")

  energy_row <- match(wellCentres, diagEnergiesData[, 1])
  E0_diag       <- diagEnergiesData[energy_row, 2]
  E1_diag       <- diagEnergiesData[energy_row, 3]
  E2_diag       <- diagEnergiesData[energy_row, 4]
  Split_diag    <- diagEnergiesData[energy_row, 5]
  Split_diag_2 <- diagEnergiesData[energy_row, 6]
}

###~~~~~~~~~~~~~~~~~###
### Acceptance rate ###
###~~~~~~~~~~~~~~~~~###

{
  acc_vec <- numeric(wellReps)

  for (i in 1:wellReps) {
    acc_vec[i] <- round(mean(accRateData[[i]]) * 100, digits = 1)
  }

  names(acc_vec) <- paste0(well_centres_vec)
  print(acc_vec)
}

###~~~~~~~~~~~~~~~~~~~~~~~~###
### Calculate ground state ### 
###~~~~~~~~~~~~~~~~~~~~~~~~###

{
  E0_vals   <- numeric(wellReps)
  E0_errors <- numeric(wellReps)

  for (i in seq_len(wellReps)) {
    
    E0_series <- E0Data[[i]]
    
    E0Split <- split(E0_series, rep(1:repeats, each = measures))
    
    E0RepeatAvg <- sapply(E0Split, mean)
    
    E0_vals[i]   <- mean(E0RepeatAvg)
    
    # Errors
    MCMC_error <- sd(E0RepeatAvg) / sqrt(length(E0RepeatAvg))

    delta_a_pc <- 0.78
    delta_beta_pc <- 1.16

    err_a    <- (delta_a_pc / 100) * mean(E0RepeatAvg)
    err_beta <- (delta_beta_pc / 100) * mean(E0RepeatAvg)
    
    # Combine errors using quadrature
    E0_errors[i] <- sqrt(MCMC_error^2 + err_a^2 + err_beta^2)
  }

  names(E0_vals)   <- well_centres_vec
  names(E0_errors) <- well_centres_vec
}

print(E0_vals, digits = 5)

###~~~~~~~~~~~~~~~~~~~~~~~###
### Excited energy states ###
###~~~~~~~~~~~~~~~~~~~~~~~###

# Log-Ratio Method

find_plateau <- function(ratios, width = 10) {
  n <- length(ratios)
  
  # Check if data is shorter than width
  if (n <= width) return(1:n)
  
  # Calculate variance for every possible window of size 'width'
  # (n - width + 1) is the number of possible windows
  variances <- sapply(1:(n - width + 1), function(i) {
    var(ratios[i:(i + width - 1)])
  })
  
  # Find the start index of the window with the minimum variance
  best_start <- which.min(variances)
  
  # Return the full range of indices for that plateau
  return(best_start:(best_start + width - 1))
}

log_ratio_energy <- function(correlator, latticeSpacing, E0, plateau_width = 10) {
  # Calculate all available log-ratios
  max_lag <- length(correlator) / 2
  raw_ratios <- numeric(max_lag)
  for (t in 1:max_lag) {
    if (correlator[t] > 0 && correlator[t + 1] > 0) {
      raw_ratios[t] <- log(correlator[t] / correlator[t + 1]) / latticeSpacing
    } else {
      raw_ratios[t] <- NA
    }
  }
  
  clean_ratios <- raw_ratios[!is.na(raw_ratios)]
  
  # Find the flattest window of the specified width
  plateau_indices <- find_plateau(clean_ratios, width = plateau_width)
  plateau_values  <- clean_ratios[plateau_indices]
  
  # Compute results
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

{
  E1_vals_log       <- numeric(wellReps)
  delta_E_vals_log  <- numeric(wellReps)
  delta_E_err_log   <- numeric(wellReps)

  noiseless_regions <- vector("list", wellReps)
  
  for (i in seq_len(wellReps)) {
    
    # Log-ratio estimation for E1
    log_ratio_results <- log_ratio_energy(
      correlator = Gx1x1Data[[i]],
      latticeSpacing = latticeSpacing,
      E0 = E0_vals[i],
      plateau_width = 20
    )
    
    E1_vals_log[i] <- log_ratio_results$E
    delta_E_vals_log[i] <- log_ratio_results$DeltaE

    names(E1_vals_log)   <- well_centres_vec
    names(delta_E_vals_log) <- well_centres_vec

    noiseless_regions[[i]] <- log_ratio_results$plateau_indices

    # Errors
    # Estimate statistical error by the standard deviation of ratios in the best plateau 
    MCMC_error  <- sd(log_ratio_results$plateau) / sqrt(length(log_ratio_results$plateau))

    delta_a_pc <- 0.78
    delta_beta_pc <- 1.16

    err_a    <- (delta_a_pc / 100) * log_ratio_results$DeltaE
    err_beta <- (delta_beta_pc / 100) * log_ratio_results$DeltaE
    
    # Combine errors using quadrature
    delta_E_err_log[i] <- sqrt(MCMC_error^2 + err_a^2 + err_beta^2)
  }
}

print(delta_E_vals_log, digits = 4)

# Exponential fit

fit_correlator <- function(correlator, latticeSpacing, noiseless_region, E0, DeltaE_guess) {
    
  correlator <- correlator[noiseless_region] # Use only noiseless region of the correlator
  
  df <- data.frame(lag = noiseless_region, correlation = correlator)
  
  fit <- nlsLM(
    correlation ~ A * exp(-DeltaE * lag * latticeSpacing) + c,
    data = df,
    start = list(A = correlator[1], c = min(correlator), DeltaE = DeltaE_guess),
    control = nls.lm.control(maxiter = 1024)
  )
  
  list(
    E = E0 + coef(fit)["DeltaE"],
    DeltaE = coef(fit)["DeltaE"],
    fit = predict(fit, newdata = df),
    df = df
  )
}

plot_correlator <- function(correlator, noiseless_region, latticeSpacing, wellCentres, exponential_fit_results = NULL) {
  n_total <- length(correlator)
  beta <- latticeSpacing * n_total
  
  full_time_axis <- (0:(n_total - 1)) * latticeSpacing
  
  df <- data.frame(
    time = full_time_axis[noiseless_region], 
    correlator = correlator[noiseless_region]
  )
  
  p <- ggplot(df, aes(x = time, y = correlator)) +
    geom_line(linewidth = 0.5) + 
    geom_point(size = 1) +
    labs(x = expression("Euclidean time" ~ tau), y = expression("G(" ~ tau ~ ")")) +
    theme_minimal(base_size = 14) +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))
  
  if (!is.null(exponential_fit_results)) {
    df_fit <- data.frame(
      time = full_time_axis[noiseless_region],
      fit_val = exponential_fit_results$fit
    )
    
    p <- p + geom_line(data = df_fit, aes(x = time, y = fit_val), 
                        color = "red", linetype = "dashed")
  }
  
  return(p)
}

{
  E1_vals_fit_naive       <- numeric(wellReps)
  delta_E_err_fit_naive   <- numeric(wellReps)
  delta_E_vals_fit_naive  <- numeric(wellReps)
  
  E1_vals_fit             <- numeric(wellReps)
  delta_E_err_fit         <- numeric(wellReps)
  delta_E_vals_fit        <- numeric(wellReps)
}

# Naively choose noiseless region

{
  noiseless_region <- 50:100

  for (i in seq_len(wellReps)) {
    exponential_fit_results <- fit_correlator(
      Gx1x1Data[[i]],
      latticeSpacing,
      noiseless_region,
      E0_vals[i],
      DeltaE_guess = Split_diag[i] 
    )

    E1_vals_fit_naive[i] <- exponential_fit_results$E
    delta_E_vals_fit_naive[i] <- exponential_fit_results$DeltaE
    delta_E_err_fit_naive[i] <- sd(exponential_fit_results$DeltaE) / sqrt(length(exponential_fit_results$DeltaE))
  }
}

delta_E_vals_fit_naive

# Choose noiseless region based on plateau from log ratio

{
  for (i in seq_len(wellReps)) {
    # Nonlinear Fit for Energy Splitting (Delta E)
    exponential_fit_results <- fit_correlator(
      Gx1x1Data[[i]],
      latticeSpacing = latticeSpacing,
      noiseless_regions[[i]],
      E0_vals[i],
      DeltaE_guess = Split_diag[i] # Could estimate splitting energy if no known values were available
    )
    
    E1_vals_fit[i] <- exponential_fit_results$E
    delta_E_vals_fit[i] <- exponential_fit_results$DeltaE
    delta_E_err_fit[i]  <- sd(exponential_fit_results$DeltaE) / sqrt(length(exponential_fit_results$DeltaE))
  }
}

delta_E_vals_fit

# Plot any of the correlator fits to check the fit
wellIndex <- 20

plot_correlator(Gx1x1Data[[wellIndex]], noiseless_regions[[wellIndex]], latticeSpacing, well_centres_vec[wellIndex],
  fit_correlator(
      Gx1x1Data[[wellIndex]],
      latticeSpacing,
      noiseless_regions[[wellIndex]],
      E0_vals[wellIndex],
      DeltaE_guess = Split_diag[wellIndex]))

if (save_png) { ggsave("Figures/DWP_exp_fit_eta_zero_point_five_window_size_twenty.png", width = 16, height = 10, dpi = 1200) }


###~~~~~~~~~~~~~~~~~~~~~~~###
### Splitting energy plot ###
###~~~~~~~~~~~~~~~~~~~~~~~###

# Create a data frame for E0
df_E0 <- data.frame(
  WellCentre = well_centres_vec,
  Energy     = E0_vals,
  Error      = E0_errors,
  Diag       = E0_diag,
  ABCs        = E0_ABCs,
  State      = "Ground State"
)

# Create a data frame for E1 
# Using the results from your Log-Ratio or Exponential Fit
df_E1 <- data.frame(
  WellCentre = well_centres_vec,
  Energy     = E1_vals_log,  # or E1_vals_fit
  Error      = sqrt(E0_errors^2 + delta_E_err_log^2),
  Diag       = E1_diag,
  ABCs        = E1_ABCs,
  State      = "1st Excited State"
)

# Combine them
plot_df_combined <- rbind(df_E0, df_E1)

ggplot(plot_df_combined, aes(x = WellCentre, y = Energy, color = State, group = State)) +
  # MCMC points and Error Bars
  geom_errorbar(aes(ymin = Energy - Error, ymax = Energy + Error), width = 0.01) +
  geom_point(size = 0.5) +
  
  # Diagonalization values as dashed lines
  geom_line(aes(y = Diag), linetype = "dashed", linewidth = 0.8, alpha = 0.5) +
  
  # ABCs values
  geom_line(aes(y = ABCs), linewidth = 0.8) +
  
  labs(
    # title = "Low-Lying Energy Spectrum of the Double Well Potential",
    # subtitle = bquote("Lattice Spacing a =" ~ .(latticeSpacing) ~ "| Inverse Temperature" ~ 
    #   beta ~ "=" ~ .(round(latticeSpacing * pathLength)) ~ "| Dashed lines = Diagonalisation | Solid lines = Semi-classical Approximation"),
    x = expression("Well Separation" ~ (eta)),
    y = "Energy (E)",
    color = "Energy Level"
  ) +
  scale_color_manual(values = c("Ground State" = "#2c3e50", "1st Excited State" = "#e74c3c")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

save_png <- TRUE

if (save_png) { ggsave("Figures/DWP_energy_spectrum.png", width = 15, height = 10, dpi = 1200) }

# Produce instanton figure for the report

# Copy paste the cout string from C++ into this vector
positions <- c(1.22811, 1.52832, 1.49021, 1.41571, 0.961768, 1.19122, 0.9546, 0.846051, 0.925141, 1.09024, 0.691262, 0.692772, 0.846204, 0.874003, 0.973807, 0.817255, 1.19464, 1.37149, 1.66603, 1.07651, 0.665956, 1.10985, 1.19542, 1.27598, 1.40074, 1.60229, 1.41072, 1.08171, 0.904115, 0.554169, 0.246528, 0.40689, 0.520851, 1.22372, 0.886372, 0.609556, 0.654539, 1.02293, 1.20372, 0.567113, 0.702199, 0.71309, 1.00189, 1.29818, 0.88723, 0.618858, 0.862358, 1.18339, 1.43116, 1.60732, 1.96511, 1.80181, 1.66373, 1.50027, 1.35927, 0.678397, 0.478808, 0.344361, 0.288119, 0.597888, 0.929044, 0.409174, 0.46435, 0.248586, 0.240017, -0.749433, -0.794577, -1.02533, -1.32361, -1.24596, -1.06513, -0.842947, -1.12487, -0.801983, -1.2046, -1.14207, -1.32984, -0.980172, -1.00385, -1.40917, -1.19175, -1.27739, -1.24336, -1.03312, -1.22493, -1.35223, -1.36932, -1.11799, -1.09601, -1.49308, -1.51708, -1.92679, -1.81992, -1.38082, -1.06052, -1.24925, -1.13222, -1.29347, -1.54426, -1.55019, -1.68911, -0.981155, -0.826883, -1.09717, -0.963735, -0.876149, -0.954911, -1.17609, -1.05865, -1.30124, -1.40704, -1.81047, -1.7651, -1.37413, -1.39044, -1.27829, -1.6689, -1.72582, -1.1289, -1.10169, -1.20334, -1.12418, -1.04066, -1.70219, -1.88639, -1.59177, -1.20228, -1.66528, -1.6055, -1.54066, -1.01199, -0.640461, -0.834568, -1.25375, -1.45562, -1.2563, -1.41776, -1.03457, -0.919233, -0.629555, -1.04076, -1.30798, -1.68327, -1.35788, -1.47487, -1.8732, -1.70886, -1.66543, -1.92222, -2.00176, -1.54993, -1.45837, -1.83612, -1.77752, -1.32215, -1.57807, -1.53441, -1.15724, -1.73176, -2.18799, -1.62336, -1.15858, -0.686137, -0.272689, -0.081936, -0.118703, -0.543992, -0.961538, -1.05379, -1.22718, -1.64245, -1.19168, -1.46999, -1.54255, -1.39312, -1.40798, -1.17618, -1.18074, -1.39286, -1.65426, -1.48059, -1.6989, -1.23399, -0.933559, -0.672166, -0.561224, -0.648794, -0.878386, -0.299315, -0.755987, -0.911297, -0.931861, -1.02026, -1.02747, -0.959071, -1.04433, -0.626058, -1.00682, -1.6379, -1.85678, -1.60383, -1.01873, -1.22281, -1.30848, -1.05561, -1.09946, -0.942813, -1.25896, -1.69456, -1.42534, -1.61618, -1.42199, -2.12937, -2.01159, -1.76654, -1.5163, -1.25219, -1.08484, -0.775807, -0.6471, 0.114379, 0.155501, 0.612678, 0.665128, 0.753684, 0.639577, 0.897472, 1.23166, 0.75444, 0.854549, 0.742459, 0.623535, 0.477715, 0.71872, 0.740067, 0.712409, 0.698456, 0.852322, 1.03546, 1.06119, 1.04849, 0.978638, 1.00106, 0.677363, 0.757804, 0.615504, 0.744954, 1.00401, 1.06618, 1.24031, 0.847568, 0.700836, 0.693529, 0.646595, 0.776235, 0.659399, 0.859068, 1.15897, 0.855017, 1.11498, 1.19787, 1.03449, 0.831773, 1.05927, 1.23089, 1.21533, 1.36913, 1.0609, 1.24488, 1.28342, 1.36737, 1.37182, 1.51418, 1.3428, 1.67839, 1.24635, 1.44194, 1.33725, 1.29749, 1.52253, 1.41181, 1.13899, 1.314, 1.12927, 0.933824, 0.91764, 1.07512, 1.36737, 1.5486, 1.69082, 1.71148, 1.79522, 1.80274, 1.09497, 1.11271, 0.999265, 1.10207, 1.28461, 1.60894, 1.54414, 1.53408, 1.29947, 1.23477, 1.26518, 1.11827, 0.987174, 1.17311, 1.05632, 0.983414, 1.34374, 1.47198, 1.37401, 1.18334, 1.20317, 1.42877, 1.32744, 1.91035, 1.55235, 1.24306, 0.858933, 0.780223, 1.26005, 1.44699, 1.38833, 1.15637, 0.555331, 0.62621, 0.792151, 0.415486, 0.629283, 0.60832, 0.902763, 1.08154, 0.801731, 0.943379, 0.604251, 1.13486, 1.499, 1.36408, 1.01323, 0.847277, 0.7362, 0.603459, 0.282474, 0.55566, 0.600056, 0.833864, 0.729355, 0.501263, 1.19513, 1.31973, 1.21708, 0.84025, 0.809774, 0.757681, 0.53272, 0.361001, 0.3696, 0.509699, 1.13068, 0.897967, 0.889386, 0.83681, 0.894613, 1.18873, 1.4772, 1.31111, 1.22256, 1.10214, 1.20643, 1.33148, 1.50443, 1.41686, 1.3373, 1.26079, 1.42095, 1.53399, 1.47653, 1.70872, 1.69623, 1.41855, 1.319, 1.34481, 1.56693, 1.1116, 1.05023, 0.99363, 1.54954, 1.56724, 1.54859, 1.60501, 1.26715, 1.45798, 1.41972, 1.48934, 1.34174, 1.40491, 1.85947, 1.67264, 1.68767, 1.42147, 1.66128, 1.66332, 1.62115, 1.44969, 1.82073, 1.5746, 0.935959, 1.29046, 1.30144, 1.07946, 0.926086, 1.24371, 0.919865, 0.959678, 0.890282, 0.494883, 0.40668, 0.269419, 0.705774, 1.03995, 1.05932, 0.369544, 0.594891, 1.33913, 1.26852, 1.08339, 1.54743, 1.23001, 1.34206, 1.40605, 1.57433, 1.56212, 1.547, 1.28023, 1.13769, 1.12138, 1.18801, 1.0472, 1.11626, 1.06166, 0.716113, 0.319511, 0.909447, 1.09824, 0.93075, 1.1055, 0.891263, 0.701796, 0.733803, 0.717993, 0.334746, 0.263842, 0.32249, 0.17851, 0.939209, 1.4337, 1.18717, 1.54779, 1.5101, 1.7799, 1.50916, 1.515, 1.15126, 1.53757, 1.61094, 1.49478, 1.62423, 1.69878, 1.68963, 1.28337, 1.13292, 1.17142, 1.00552, 0.7865, 1.00752, 1.3144, 1.12108, 1.42405, 1.46268, 1.25337, 1.33368, 1.0103, 0.826797, 0.555796, 0.842378, 0.87257, 1.08469, 0.746488, 1.08138, 1.00486, 1.00237, 1.10258, 0.990965, 1.18924, 0.904281, 0.787234, 0.731169, 1.17486, 1.16662)

a <- 0.075
path_length <- 500
tau <- seq(0, (path_length - 1) * a, by = a)

# Define physical parameters from simulation
wellCentres <- 1.4
lambda <- 1.0
omegaDWP <- wellCentres * sqrt(2 * lambda)

# Define as many times as there are instantons
tau_0 <- 4.9
tau_1 <- 16.7

# Simply multiply this by further tanh functions with tau_x as new instantons/anti-instantons
instanton_sol <- wellCentres * tanh(omegaDWP * (tau - tau_0)) * tanh(omegaDWP * (tau - tau_1))

df <- data.frame(
  tau = tau,
  positions = positions,
  instanton = instanton_sol
)

ggplot(df, aes(x = tau)) +
  geom_line(aes(y = positions), color = "black", linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = instanton), color = "red", linewidth = 1, linetype = "dashed") +
  labs(x = expression("Euclidean time" ~ tau), y = "Position x") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), legend.position = "none")

if (save_png) { ggsave("Figures/anti-instanton-instanton.png", width = 15, height = 10, dpi = 1200) }

# nolint end