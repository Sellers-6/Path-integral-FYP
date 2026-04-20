# nolint start

#####################################################################################################################################
#################################################### Double Well Potential ##########################################################
#####################################################################################################################################

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Libraries and helper functions ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

{
  library(ggplot2) 
  library(gridExtra)    # Side by side plots
  library(dplyr)
  library(rhdf5)        # Read h5 files
  library(minpack.lm)   # Better exponential fitting
  library(expm)         # Matrix exponentials

  fit_correlator_exp <- function(corr, noiseless_region, A_guess, DeltaE_guess) {
    
    lag <- 0:(noiseless_region - 1)
    corr_slice <- corr[1:noiseless_region]
    
    dfCorr <- data.frame(
      lag = lag,
      correlation = corr_slice
    )
    
    c_guess <- min(corr_slice, na.rm = TRUE)
    
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
        geom_point(aes(y = correlation, color = "MCMC Correlator"), size = 0.5) +
        geom_line(aes(y = fit, color = "Exponential fit")) +
        labs(
          title = "Correlator with Exponential Fit",
          y = "G(t)",
          x = "t"
        ) +
        scale_color_manual(
          name = "",
          values = c("MCMC Correlator" = "black", "Exponential fit" = "red")
        )
    
    list(
      coefficients = coef(fit),
      fit_object = fit,
      data = dfCorr,
      plot = p
    )
  }

  log_ratio_energy <- function(corr, latticeSpacing, noiseless_region, E0) {
    
    if (length(corr) < noiseless_region + 1) {
      stop("Correlator too short for chosen noiseless_region")
    }
    
    RHS <- noiseless_region
    LHS <- 1
    
    ratios <- numeric(RHS - LHS + 1)
    valid_count <- 0
    
    # Compute log ratios
    for (i in LHS:RHS) {
      
      if (corr[i] > 0 && corr[i + 1] > 0) {
        
        r <- log(corr[i] / corr[i + 1]) / latticeSpacing
        
        ratios[i - LHS + 1] <- r
        valid_count <- valid_count + 1
        
      } else {
        ratios[i - LHS + 1] <- NA
      }
    }
    
    # Average over valid log ratios to get DeltaE
    valid_ratios <- ratios[!is.na(ratios)]
    
    if (length(valid_ratios) == 0) {
      stop("No valid log-ratio points found")
    }
    
    DeltaE <- mean(valid_ratios)
    excited_energy <- E0 + DeltaE
    
    df <- data.frame(
      lag = LHS:RHS,
      ratio = ratios
    )
    
    p <- ggplot(df, aes(x = lag * latticeSpacing, y = ratio)) +
      geom_line(color = "blue") +
      geom_hline(
        yintercept = DeltaE,
        linetype = "dashed",
        color = "red"
      ) +
      labs(
        title = "Log-Ratio Splitting Energy Estimate",
        x = "t",
        y = "Splitting Energy"
      )
    
    list(
      excited_energy = excited_energy,
      DeltaE = DeltaE,
      ratios = ratios,
      mean_ratios = valid_ratios,
      plot = p
    )
  }
}

# Overwrite figures?
save_png <- FALSE

###~~~~~~~~~~~###
### Read data ###
###~~~~~~~~~~~###

{
  dataFile <- "data2.h5"
  
  E0ThermData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/E0Therm"))))
  accRateThermData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/accRateTherm"))))
  E0Data <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/E0"))))
  accRateData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/accRate"))))
  histogramData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/histogram"))))
  instantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/instantons"))))
  antiInstantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/antiInstantons"))))
  Gx1x1Data <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/Gx1x1"))))
  Gx2x2Data <- as.numeric(unlist(h5read(dataFile, paste0("/DWP/Gx2x2"))))
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
    energy_row    <- as.integer(round((wellCentres - 0.95) * 20))
    E0_diag       <- diagEnergiesData[energy_row, 2]
    E1_diag       <- diagEnergiesData[energy_row, 3]
    E2_diag       <- diagEnergiesData[energy_row, 4]
    Split_diag    <- diagEnergiesData[energy_row, 5]
    Split_diag_2 <- diagEnergiesData[energy_row, 6]

    # WKB analysis
    omegaDWP  <- sqrt(8 * lambdaDWP * wellCentres^2)
    S_inst    <- (2/3) * omegaDWP * (wellCentres^2) 
    alpha     <- 1 / 12 # A complicated calculation performed in Zinn-Justin 1993 (or ABCs of Instantons) gives this value
    K         <- omegaDWP * sqrt(S_inst / (2 * pi)) * (alpha ^ -0.5) # A prefactor for the splitting energy
    Split_WKB <- K * exp(-S_inst)
    E0_WKB    <- 0.5 * omegaDWP - (Split_WKB / 2)
    E1_WKB    <- 0.5 * omegaDWP + (Split_WKB / 2)
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

  E0_error <- abs(((E0 - E0_diag) / E0_diag) * 100)
}

E0; E0_diag; mean(E0RepeatAvg) + E0StandardError; mean(E0RepeatAvg) - E0StandardError ; E0_error

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
      labs(title = "Ground State Energy Histogram for the DWP",
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
ggplot(wave_df, aes(x = x)) +
    # MCMC wave function
    geom_point(aes(y = psi, color = "MCMC"), size = 1) +

    # Diagonalisation wave function
    geom_line(data = diagWFData, aes(x = x, y = psi0, color = "Diagonalised"), linewidth = 0.7) +

    # WKB approximation
    geom_line(aes(y = psiAnalytical, color = "WKB"), linetype = "dashed") +
    
    labs(title = "Double-Well Potential Ground State Wave Function", x = "Position x") +

    # coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 0.8)) + 

    scale_y_continuous(name = expression("Wavefunction " ~ psi[0](x))) +

    scale_color_manual(name = "", values = c("MCMC" = "red", "Diagonalised" = "black", "WKB" = "blue")) + 

    theme_minimal(base_size = 14)
   
if (save_png) { ggsave("Figures/DWP_wave_functions.png", width = 8, height = 4, dpi = 1200) }

###~~~~~~~~~~~~~~~~~~~~~~~###
### Excited energy states ###
###~~~~~~~~~~~~~~~~~~~~~~~###

### Exponential fit ###

{
  Gx1x1Fit <- fit_correlator_exp(corr = Gx1x1Data, noiseless_region = 50, A_guess = E0_diag, DeltaE_guess = Split_diag)
  Gx2x2Fit <- fit_correlator_exp(corr = Gx2x2Data, noiseless_region = 30, A_guess = E0_diag, DeltaE_guess = Split_diag_2)
}

Gx1x1Fit$coefficients
print(Gx1x1Fit$plot) 

if (save_png) { ggsave("Figures/DWP_x1x1_correlator.png", width = 7, height = 4, dpi = 1200) }

Gx2x2Fit$coefficients
print(Gx2x2Fit$plot) 
 
if (save_png) { ggsave("Figures/DWP_x2x2_correlator.png", width = 7, height = 4, dpi = 1200) }

### Log ratio method ### 

{
  Gx1x1_log_ratios <- log_ratio_energy(corr = Gx1x1Data, latticeSpacing = latticeSpacing, noiseless_region = 500, E0 = E0)
  Gx2x2_log_ratios <- log_ratio_energy(corr = Gx2x2Data, latticeSpacing = latticeSpacing, noiseless_region = 50, E0 = E0)
}

Gx1x1_log_ratios$excited_energy
Gx1x1_log_ratios$DeltaE
print(Gx1x1_log_ratios$plot)

Gx2x2_log_ratios$excited_energy
Gx2x2_log_ratios$DeltaE
print(Gx2x2_log_ratios$plot)

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

# Comparing energies against diagonalisation method and WKB

E0; E0_diag; E0_WKB

E1; E1_diag; E1_WKB

E1 - E0; Split_diag; Split_WKB 

#####################################################################################################################################
################################################# Vary Well Centres and Beta ########################################################
#####################################################################################################################################

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Libraries and helper functions ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

{
  library(ggplot2) 
  library(gridExtra)    # Side by side plots
  library(dplyr)
  library(rhdf5)        # Read h5 files
  library(minpack.lm)   # Better exponential fitting
  library(expm)         # Matrix exponentials

  read_h5 <- function(name) as.numeric(unlist(h5read(dataFile, paste0(base, name))))

  fit_correlator <- function(correlator, latticeSpacing, noiseless_region, E0, DeltaE_guess) {
    
    lag <- 0:(noiseless_region - 1)
    correlator <- correlator[1:noiseless_region] # Use only noiseless region of the correlator
    
    df <- data.frame(lag = lag, correlation = correlator)
    
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

  log_ratio_energy <- function(correlator, latticeSpacing, E0, max_lag) {
    
    ratios <- numeric(max_lag)
    
    for (i in 1:max_lag) {
      if (correlator[i] > 0 && correlator[i + 1] > 0) { 
        ratios[i] <- log(correlator[i] / correlator[i + 1]) / latticeSpacing 
      } 
      else { ratios[i] <- NA }
    }
    
    ratios <- ratios[!is.na(ratios)] # Reject NA log results
    
    DeltaE <- mean(ratios)
    E <- E0 + DeltaE
    
    list(
      E = E,
      DeltaE = DeltaE,
      ratios = ratios
    )
  }

  plot_correlator <- function(correlator, latticeSpacing, exponential_fit_results = NULL) {
    n <- length(correlator)
    df <- data.frame(time = 0:(n - 1) * latticeSpacing, correlator = correlator)
    
    p <- ggplot(df, aes(x = time, y = correlator)) +
      geom_line(size = 0.5) +
      geom_point(size = 1) +
      labs(x = "Euclidean time t", y = "G(t)", title = "DWP correlator") +
      theme_minimal(base_size = 14)
    
    if (!is.null(exponential_fit_results)) {
      df_fit <- data.frame(time = (0:(length(exponential_fit_results$fit) - 1)) * latticeSpacing,
                          fit_val = exponential_fit_results$fit)
      
      p <- p + geom_line(data = df_fit, aes(x = time, y = fit_val), color = "red", linetype = "dashed")
    }
    
    return(p)
  }
}

# Overwrite figures?
save_png <- FALSE

###~~~~~~~~~~~###
### Read data ###
###~~~~~~~~~~~###

{
  dataFile <- "data2.h5"
  
  # well_centres_vec <- c(1.0, 1.1, 1.2, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.6)
  # betas_vec <- c(1000, 750, 500, 250, 100, 75, 50, 25, 10)
  well_centres_vec <- c(1.3, 1.4, 1.5, 1.6)
  betas_vec <- c(1000, 500, 100, 50, 10)
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
  E0Data             <- vector("list", wellReps)
  accRateData        <- vector("list", wellReps)
  histogramData      <- vector("list", wellReps)
  Gx1x1Data          <- vector("list", wellReps)
  Gx2x2Data          <- vector("list", wellReps)
  instantonsData     <- vector("list", wellReps)
  antiInstantonsData <- vector("list", wellReps)
  headerData         <- vector("list", wellReps)
  
  # Read data
  for (i in seq_along(well_centres_vec)) {
    
    # Initialize sublists for each well center
    E0Data[[i]]             <- vector("list", betaReps)
    accRateData[[i]]        <- vector("list", betaReps)
    histogramData[[i]]      <- vector("list", betaReps)
    instantonsData[[i]]     <- vector("list", betaReps)
    antiInstantonsData[[i]] <- vector("list", betaReps)
    Gx1x1Data[[i]]          <- vector("list", betaReps)
    Gx2x2Data[[i]]          <- vector("list", betaReps)
    headerData[[i]]         <- vector("list", betaReps)
    
    for (j in seq_along(betas_vec)) {
      base <- full_paths[i, j]
      
      # Read everything into the nested list structure
      E0Data[[i]][[j]]             <- read_h5("E0")
      accRateData[[i]][[j]]        <- read_h5("accRate")
      histogramData[[i]][[j]]      <- read_h5("histogram")
      instantonsData[[i]][[j]]     <- read_h5("instantons")
      antiInstantonsData[[i]][[j]] <- read_h5("antiInstantons")
      Gx1x1Data[[i]][[j]]          <- read_h5("Gx1x1")
      Gx2x2Data[[i]][[j]]          <- read_h5("Gx2x2")
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

  Split_WKB <- K * exp(-S_inst)
  E0_WKB <- 0.5 * omegaDWP - (Split_WKB / 2)
  E1_WKB <- 0.5 * omegaDWP + (Split_WKB / 2)

  # Diagonalisation data
  diagEnergiesData <- read.csv("DWP diagonalisation/DWP diagonalisation/energies.csv")

  energy_row <- as.integer(round((wellCentres - 0.95) * 20))

  E0_diag       <- diagEnergiesData[energy_row, 2]
  E1_diag       <- diagEnergiesData[energy_row, 3]
  E1_diag       <- diagEnergiesData[energy_row, 4]
  Split_diag    <- diagEnergiesData[energy_row, 5]
  Split_diag_2 <- diagEnergiesData[energy_row, 6]
}

###~~~~~~~~~~~~~~~~~###
### Acceptance rate ###
###~~~~~~~~~~~~~~~~~###

acc_matrix <- matrix(0, nrow = wellReps, ncol = betaReps)

for (i in 1:wellReps) {
  for (j in 1:betaReps) {
    acc_matrix[i, j] <- mean(accRateData[[i]][[j]]) * 100
  }
}

rownames(acc_matrix) <- paste0("Well Centre: ", well_centres_vec)
colnames(acc_matrix) <- paste0("Beta: ", betas_vec)

print(acc_matrix)

###~~~~~~~~~~~~~~~~~~~~~~~~###
### Calculate ground state ### 
###~~~~~~~~~~~~~~~~~~~~~~~~###

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

# Ground state energies
E0_vals[,]
# Error bars
E0_errors[,]

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

ggplot(plot_df, aes(x = Beta, y = Energy, color = Well_Color)) +
  # Diagonalised values
  geom_hline(aes(yintercept = E0_diag, color = Well_Color), linetype = "dashed") +
  
  # MCMC ground state energy
  geom_line(aes(group = WellCentre), size = 0.7) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Energy - Error, ymax = Energy + Error), width = 0.05) +
  
  # Apply log scale to X-axis as Beta spans orders of magnitude
  scale_x_log10(breaks = betas_vec) +
  
  # Labels and theme
  labs(title = "Ground State Energy Convergence", x = "Inverse Temperature Beta", 
      y = expression("Ground State Energy " ~ E[0]), color = "Well Centre") +
  theme_minimal(base_size = 14)

###~~~~~~~~~~~~~~~~~~~~~~~###
### Excited energy states ###
###~~~~~~~~~~~~~~~~~~~~~~~###

noiseless_region <- 20

# Exponential fit method

E1_vals_fit       <- matrix(NA, wellReps, betaReps)
delta_E_vals_fit  <- matrix(NA, wellReps, betaReps)
delta_E_err_fit  <- matrix(NA, wellReps, betaReps)

for (i in seq_len(wellReps)) {
  for (j in seq_len(betaReps)) {

    # Nonlinear Fit for Energy Splitting (Delta E)
    exponential_fit_results <- fit_correlator(
      Gx1x1Data[[i]][[j]],
      latticeSpacing,
      noiseless_region,
      E0_vals[i, j],
      DeltaE_guess = Split_diag[i] # Could estimate splitting energy if no known values were available
    )
    
    E1_vals_fit[i, j] <- exponential_fit_results$E
    delta_E_vals_fit[i, j] <- exponential_fit_results$DeltaE
    delta_E_err_fit[i, j]  <- sd(exponential_fit_results$DeltaE) / sqrt(length(exponential_fit_results$DeltaE))
  }
}

E1_vals_fit
delta_E_vals_fit

# Log ratio method

E1_vals_log       <- matrix(NA, wellReps, betaReps)
delta_E_vals_log  <- matrix(NA, wellReps, betaReps)
delta_E_err_log  <- matrix(NA, wellReps, betaReps)

for (i in seq_len(wellReps)) {
  for (j in seq_len(betaReps)) {
    # Log-ratio estimation for E1
    log_ratio_results <- log_ratio_energy(
      corr = Gx1x1Data[[i]][[j]],
      latticeSpacing = latticeSpacing,
      E0 = E0_vals[i, j],
      max_lag = noiseless_region
    )
    
    E1_vals_log[i, j] <- log_ratio_results$E
    delta_E_vals_log[i, j] <- log_ratio_results$DeltaE
    delta_E_err_log[i, j]  <- sd(log_ratio_results$DeltaE) / sqrt(length(log_ratio_results$DeltaE))
  }
}

E1_vals_log
delta_E_vals_log

wellIndex <- 1
betaIndex <- 1

plot_correlator(Gx1x1Data[[wellIndex]][[betaIndex]][1:noiseless_region], latticeSpacing, 
  fit_correlator(
      Gx1x1Data[[wellIndex]][[betaIndex]],
      latticeSpacing,
      noiseless_region,
      E0_vals[wellIndex, betaIndex],
      DeltaE_guess = Split_diag[betaIndex]))

plot_fit_df <- expand.grid(
  WellCentre = well_centres_vec,
  Beta = betas_vec
)

# Flatten matrices into the data frame
plot_fit_df$E1        <- as.vector(E1_vals)
plot_fit_df$E1_Err    <- as.vector(E1_err)
plot_fit_df$Splitting <- as.vector(Split_fit)

# 2. Process with dplyr
plot_fit_df <- plot_fit_df %>%
  mutate(
    Well_Factor = as.factor(WellCentre)
  )

# 3. Create the Plot
ggplot(plot_fit_df, aes(x = Beta, y = Splitting, color = Well_Factor)) +
  # Add the simulation lines and points
  geom_line(aes(group = WellCentre), size = 0.8) +
  geom_point(size = 2.5) +
  
  labs(
    title = "Energy Splitting Convergence with increasing Beta",
    x = expression("Inverse Temperature" ~ Beta),
    y = expression(Delta * E),
    color = "Well Centre") +
  theme_minimal(base_size = 14)


# # Exponential fit
# E1_vals <- matrix(NA, wellReps, latReps)
# E1_err  <- matrix(NA, wellReps, latReps)

# noiseless_region <- 10

# for (i in seq_len(wellReps)) {
#   for (j in seq_len(latReps)) {
    
#     latticeSpacing <- lattice_spacings_vec[j]
#     E0_here <- E0_vals[i, j]
    
#     # G_1(n) = G(x(t),x(t + n)) correlator gives E1
#     corr1 <- Gx1x1Data[i, j, ]
    
#     res1 <- log_ratio_energy(
#       corr = corr1,
#       latticeSpacing = latticeSpacing,
#       E0 = E0_here,
#       max_lag = noiseless_region
#     )
    
#     E1_vals[i, j] <- res1$E
    
#     # crude error estimate
#     E1_err[i, j] <- sd(res1$ratios) / sqrt(length(res1$ratios))
#   }
# }

# Split_fit <- matrix(NA, wellReps, latReps)

# for (i in seq_len(wellReps)) {
#   for (j in seq_len(latReps)) {
    
#     latticeSpacing <- lattice_spacings_vec[j]
#     corr <- Gx1x1Data[i, j, ]
    
#     fit_res <- fit_correlator(
#       corr,
#       latticeSpacing,
#       noiseless_region,
#       DeltaE_guess = 0.1 # Put the diagonalisation splitting energy here as a guess
#     )
    
#     Split_fit[i, j] <- fit_res$DeltaE
#   }
# }

# plot_correlator <- function(corr, latticeSpacing, fit_res = NULL) {
  
#   n <- length(corr)
#   df <- data.frame(
#     lag = 0:(n - 1),
#     corr = corr
#   )
  
#   p <- ggplot(df, aes(x = lag * latticeSpacing, y = corr)) +
#     geom_line(size = 0.5) +
#     labs(x = "t", y = "G(t)")
  
#   if (!is.null(fit_res)) {
#     df$fit <- fit_res$fit
#     p <- p + geom_line(aes(y = fit), color = "red")
#   }
  
#   p
# }

# plot_correlator(Gx1x1Data[10,1,][1:10], 0.1)

# E1_vals[,1]

# E2_vals[,1]

# Split_fit[,1]


# nolint end