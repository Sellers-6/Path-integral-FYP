# nolint start

library(ggplot2)
library(gridExtra)  # side by side plots
library(dplyr)


library(rhdf5)

h5ls("simulations.h5")

vec <- h5read("simulations.h5", "/E0Therm/Periodic/QHO")
head(vec)


### Reading files ###

# Define boundary conditions and system types
BCs <- c("Periodic", "Dirichlet") 
systems <- c("QHO", "DWP") 

# Create empty lists to store data
therm <- list() 
E0evolution <- list()
waveFunction <- list()
correlation <- list()
E0 <- list()
E1 <- list()
accRate <- list()



# Helper function to read a CSV safely
readCsvSafe <- function(fileName, header = TRUE) {
  if (file.exists(fileName)) {
    read.csv(fileName, header = header)
  } else {
    warning(paste("File not found:", fileName))
    return(NULL)
  }
}

# Loop over all combinations
for (bc in BCs) {
  for (sys in systems) {
    name <- paste0(bc, sys)
    
    thermFile <- paste0("csv/E0Thermalisation", bc, sys, ".csv")
    E0evolFile <- paste0("csv/E0Evolution", bc, sys, ".csv")
    waveFile <- paste0("csv/waveFunction", bc, sys, ".csv")
    corrFile <- paste0("csv/correlation", bc, sys, ".csv")
    E0File <- paste0("csv/E0", bc, sys, ".csv")
    E1File <- paste0("csv/E1", bc, sys, ".csv")
    accRateFile <- paste0("csv/accRate", bc, sys, ".csv")
     


    # Read each CSV safely and store directly in the lists
    therm[[name]] <- readCsvSafe(thermFile, header = TRUE)
    E0evolution[[name]] <- readCsvSafe(E0evolFile, header = TRUE)
    waveFunction[[name]] <- readCsvSafe(waveFile, header = FALSE)
    correlation[[name]] <- readCsvSafe(corrFile, header = TRUE)
    E0[[name]] <- readCsvSafe(E0File, header = TRUE)
    E1[[name]] <- readCsvSafe(E1File, header = TRUE)
    accRate[[name]] <- readCsvSafe(accRateFile, header = TRUE)

  }
}

# bc <- "Periodic"    # or "Dirichlet"
bc <- "Dirichlet"   # or "Periodic"
sys <- "QHO"        # or "DWP"
# sys <- "DWP"        # or "QHO"

a <- 2              # DWP variables
lambda <- 1/12

name <- paste0(bc, sys)      # e.g., "PeriodicQHO"

# Thermalisation

ggplot(therm[[name]], aes(x = as.numeric(row.names(therm[[name]])), y = E0)) +
    geom_point() +
    labs(
      x = "Measurement index",
      y = "E0",
      title = paste("Thermalisation of", bc, sys)
    )

# Decorrelation

ggplot(E0evolution[[name]], aes(x = as.numeric(row.names(E0evolution[[name]])), y = E0)) +
  geom_point() + 
  geom_line() + 
  labs(
      x = "Measurement index",
      y = "E0",
      title = paste("Decorrelation of", bc, sys)
    ) 

# Wave function

matrixOfPositions <- as.numeric(as.matrix(waveFunction[[name]])) 
xVec <- as.vector(matrixOfPositions) 

bins <- 100

hist(xVec, breaks = bins, main = "Histogram of positions",
     xlab = "position", ylab = "count")
h <- hist(xVec, breaks = 100, plot = FALSE)


measures <- sum(h$counts)                             # Total counts which equals measures * lattice size
positionRange <- range(h$mids)[2] - range(h$mids)[1]  # Range of positions 
binWidth <- positionRange / bins                      # Width of each bin 
normFactor <- measures * binWidth                     # Normalization factor 
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

ggplot(correlation[[name]], aes(x = as.numeric(row.names(correlation[[name]])), y = Correlation)) + 
  geom_point() +
  labs(title = "Correlation Function", x = "Position", y = "Correlation")

L <- 5000 
t <- seq(0, L, length.out = L)

correlator <- function(t) {
  0.5 * exp(-t) + 0.5 * exp(t - (L - 1))
}

gAnalytic <- data.frame( 
  t = t,
  G = sapply(0:(L - 1), correlator)
)

ggplot(correlation[[name]], aes(x = as.numeric(row.names(correlation[[name]])), y = Correlation)) + 
  geom_point() +
  geom_line(data = gAnalytic, aes(x = t, y = G), color = "red") +
  labs(x = "Distance", y = "Correlation") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))

# Ground state energy histogram

e0Vec <- as.numeric(E0[[name]]$E0)

bins <- 10 

hE0 <- hist(e0Vec, breaks = bins, plot = FALSE) # Compute histogram breaks and bin width
binWidth <- diff(hE0$breaks)[1]

E0Norm <- hE0$counts / (sum(hE0$counts) * binWidth) # Normalised histogram counts (probability density)

continuousE0 <- seq(min(e0Vec), max(e0Vec), length.out = 200) # Continuous x values for normal curve
normDist <- dnorm(continuousE0, mean = mean(e0Vec), sd = sd(e0Vec))
dx <- diff(continuousE0)[1]
normDist <- normDist / sum(normDist * dx)  # Normalisation

# Histogram and normal curve
histPlot <- ggplot() +
  geom_histogram(aes(x = e0Vec, y = after_stat(density)),
                 bins = bins, fill = "skyblue", color = "black") +
  geom_line(aes(x = continuousE0, y = normDist),
            color = "red", linewidth = 1) +
  labs(title = "Ground State Energy Histogram",
       x = "Energy", y = "Probability Density") +
  theme_minimal()

# QQ plot
shapiroTest <- shapiro.test(e0Vec)
qqPlot <- ggplot(data.frame(e0Vec), aes(sample = e0Vec)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = paste0("QQ Plot (Shapiro-Wilk p = ", round(shapiroTest$p.value, 4), ")"),
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

# Combined plots side by side
# grid.arrange(histPlot, qqPlot, ncol = 2)
histPlot


# First excited energy state histogram

e1Vec <- as.numeric(E1[[name]]$E1)

bins <- 5

hist(e1Vec, breaks = bins, main = "Histogram of first excited state energies",
     xlab = "Energy", ylab = "count")


# nolint end
