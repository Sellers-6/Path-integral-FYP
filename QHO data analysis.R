# nolint start

#####################################################################################################################################
################################################ Quantum Harmonic Oscillator ########################################################
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
        c = c_guess
      ),
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
  
  E0ThermData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/E0Therm"))))
  accRateThermData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/accRateTherm"))))
  E0Data <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/E0"))))
  accRateData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/accRate"))))
  histogramData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/histogram"))))
  instantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/instantons"))))
  antiInstantonsData <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/antiInstantons"))))
  Gx1x1Data <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/Gx1x1"))))
  Gx2x2Data <- as.numeric(unlist(h5read(dataFile, paste0("/QHO/Gx2x2"))))
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

  E0_error <- abs(((E0 - 0.5) / 0.5) * 100) # Percent error
}

E0; mean(E0RepeatAvg) + E0StandardError; mean(E0RepeatAvg) - E0StandardError ; E0_error

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
ggplot(wave_df, aes(x = x)) +
    # MCMC wave function
    geom_point(aes(y = psi, color = "MCMC"), size = 1) +

    # Exact wave function
    geom_line(aes(y = psiAnalytical, color = "Exact")) +

    labs(title = paste("Quantum Harmonic Oscillator Ground State Wave Function"), x = "Position", y = "Psi") +

    scale_y_continuous(name = expression("Wavefunction " ~ psi[0](x))) +

    scale_color_manual(name = "", values = c("MCMC" = "red", "Exact" = "black")) + 

    theme_minimal(base_size = 14)
   
if (save_png) { ggsave("Figures/QHO_wave_functions.png", width = 8, height = 4, dpi = 1200) }

###~~~~~~~~~~~~~~~~~~~~~~~###
### Excited energy states ###
###~~~~~~~~~~~~~~~~~~~~~~~###

### Exponential fit ###

{
  Gx1x1Fit <- fit_correlator_exp(corr = Gx1x1Data, noiseless_region = 50, A_guess = 0.5, DeltaE_guess = 1.0)
  Gx2x2Fit <- fit_correlator_exp(corr = Gx2x2Data, noiseless_region = 30, A_guess = 0.5, DeltaE_guess = 2.0)
}

Gx1x1Fit$coefficients
print(Gx1x1Fit$plot) 

if (save_png) { ggsave("Figures/QHO_x1x1_correlator.png", width = 7, height = 4, dpi = 1200) }

Gx2x2Fit$coefficients
print(Gx2x2Fit$plot) 
 
if (save_png) { ggsave("Figures/QHO_x2x2_correlator.png", width = 7, height = 4, dpi = 1200) }

### Log ratio method ### 

{
  Gx1x1_log_ratios <- log_ratio_energy(corr = Gx1x1Data, latticeSpacing = latticeSpacing, noiseless_region = 25, E0 = E0)
  Gx2x2_log_ratios <- log_ratio_energy(corr = Gx2x2Data, latticeSpacing = latticeSpacing, noiseless_region = 10, E0 = E0)
}

Gx1x1_log_ratios$excited_energy
print(Gx1x1_log_ratios$plot)

Gx2x2_log_ratios$excited_energy
print(Gx2x2_log_ratios$plot)

# nolint end