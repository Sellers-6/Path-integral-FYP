library(ggplot2)
library(dplyr)

### Periodic boundary conditions ###

## Quantum harmonic oscillator ##

thermPQHO <- read.csv("csv/E0ThermalisationPeriodicQHO.csv") #nolint

ggplot(thermPQHO, aes(x = as.numeric(row.names(thermPQHO)), y = E0)) +
  geom_point()


E0evolutionPQHO <- read.csv("csv/E0EvolutionPeriodicQHO.csv") #nolint

ggplot(E0evolutionPQHO, aes(x = as.numeric(row.names(E0evolutionPQHO)), y = E0)) + #nolint
  geom_point() + geom_line()

waveFunctionPQHO <- read.csv("csv/waveFunctionPeriodicQHO.csv", header = FALSE) #nolint

matrixOfPositions <- as.matrix(waveFunctionPQHO) #nolint

xVec <- as.vector(matrixOfPositions) #nolint

bins <- 50

hist(xVec, breaks = bins, main = "Histogram of positions",
     xlab = "position", ylab = "count")
h <- hist(xVec, breaks = 200, plot = FALSE)

measures <- sum(h$counts)  # total counts which equals measures * lattice size
positionRange <- range(h$mids)[2] - range(h$mids)[1]  # range of positions # nolint
binWidth <- positionRange / bins  # width of each bin # nolint
normFactor <- measures * binWidth  # normalization factor # nolint
psi <- sqrt(h$counts / normFactor)  # normalized wave function # nolint

sum(psi ^ 2) * binWidth  # should be approximately 1 # nolint

psiAnalytical <- exp(-(h$mids ^ 2) / 2) # Obtaining the analytical wavefunction #nolint

psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) # Normalising the analytical wavefunction #nolint

sum(psiAnalytical ^ 2) * binWidth

ggplot(data.frame(x = h$mids, psi = psi), aes(x = x, y = psi)) +
  geom_line() +
  geom_point() +
  geom_line(aes(x = h$mids, y = psiAnalytical)) +
  labs(title = "Wave Function", x = "Position", y = "Psi")

correlationPQHO <- read.csv("csv/correlationPeriodicQHO.csv") #nolint

head(correlationPQHO)
ggplot(correlationPQHO, aes(x = as.numeric(row.names(correlationPQHO)), y = Correlation)) + #nolint
  geom_point()

min(correlationPQHO$Correlation)
which(correlationPQHO$Correlation == min(correlationPQHO$Correlation))
L <- 5000 #nolint

t <- seq(0, L, length.out = L)

correlator <- function(t) {
  0.5 * exp(-t) + 0.5 * exp(t - (L - 1))
}

gAnalytic <- data.frame( #nolint
  t = t,
  G = sapply(0:(L - 1), correlator)
)

ggplot(correlationPQHO, aes(x = as.numeric(row.names(correlationPQHO)), y = Correlation)) + #nolint
  geom_point() +
  geom_line(data = gAnalytic, aes(x = t, y = G), color = "red") +
  labs(x = "Distance", y = "Correlation") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))

## Double Well Potential ##

thermPDWP <- read.csv("csv/E0ThermalisationPeriodicDWP.csv") #nolint

ggplot(thermPDWP, aes(x = as.numeric(row.names(thermPDWP)), y = E0)) +
  geom_point()


E0evolutionPDWP <- read.csv("csv/E0EvolutionPeriodicDWP.csv") #nolint
ggplot(E0evolutionPDWP, aes(x = as.numeric(row.names(E0evolutionPDWP)), y = E0)) + #nolint
  geom_point() + geom_line()

waveFunctionPDWP <- read.csv("csv/waveFunctionPeriodicDWP.csv") #nolint
x <- waveFunctionPDWP$Position # Temporary variable to store the positions #nolint

bins <- 50

hist(x, breaks = bins, main = "Histogram of positions",
     xlab = "position", ylab = "count")
h <- hist(x, breaks = 200, plot = FALSE)

measures <- sum(h$counts)  # total counts which equals measures * lattice size
positionRange <- range(h$mids)[2] - range(h$mids)[1]  # range of positions # nolint
binWidth <- positionRange / bins  # width of each bin # nolint
normFactor <- measures * binWidth  # normalization factor # nolint
psi <- sqrt(h$counts / normFactor)  # normalized wave function # nolint

sum(psi ^ 2) * binWidth  # should be approximately 1 # nolint

psiAnalytical <- exp(-(h$mids ^ 2) / 2) # Obtaining the analytical wavefunction #nolint

psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) # Normalising the analytical wavefunction #nolint

sum(psiAnalytical ^ 2) * binWidth

ggplot(data.frame(x = h$mids, psi = psi), aes(x = x, y = psi)) +
  geom_line() +
  geom_point() +
  geom_line(aes(x = h$mids, y = psiAnalytical)) +
  labs(title = "Wave Function", x = "Position", y = "Psi")

correlationPDWP <- read.csv("csv/correlationPeriodicDWP.csv") #nolint

head(correlationPDWP)
ggplot(correlationPDWP, aes(x = as.numeric(row.names(correlationPDWP)), y = Correlation)) + #nolint
  geom_point()

min(correlationPDWP$Correlation)
which(correlationPDWP$Correlation == min(correlationPDWP$Correlation))
L <- 5000 #nolint

t <- seq(0, L, length.out = L)

correlator <- function(t) {
  0.5 * exp(-t) + 0.5 * exp(t - (L - 1))
}

gAnalytic <- data.frame( #nolint
  t = t,
  G = sapply(0:(L - 1), correlator)
)

ggplot(correlationPDWP, aes(x = as.numeric(row.names(correlationPDWP)), y = Correlation)) + #nolint
  geom_point() +
  geom_line(data = gAnalytic, aes(x = t, y = G), color = "red") +
  labs(x = "Distance", y = "Correlation") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))



### Dirichlet boundary conditions ###


## Quantum harmonic oscillator ##

thermDQHO <- read.csv("csv/E0ThermalisationDirichletQHO.csv") #nolint

ggplot(thermDQHO, aes(x = as.numeric(row.names(thermDQHO)), y = E0)) +
  geom_point()


E0evolutionDQHO <- read.csv("csv/E0EvolutionDirichletQHO.csv") #nolint
ggplot(E0evolutionDQHO, aes(x = as.numeric(row.names(E0evolutionDQHO)), y = E0)) + #nolint
  geom_point() + geom_line()

waveFunctionDQHO <- read.csv("csv/waveFunctionDirichletQHO.csv") #nolint

x <- waveFunctionDQHO$Position

bins <- 50

hist(x, breaks = bins, main = "Histogram of positions",
     xlab = "position", ylab = "count")
h <- hist(x, breaks = 200, plot = FALSE)

measures <- sum(h$counts)  # total counts which equals measures * lattice size
positionRange <- range(h$mids)[2] - range(h$mids)[1]  # range of positions # nolint
binWidth <- positionRange / bins  # width of each bin # nolint
normFactor <- measures * binWidth  # normalization factor # nolint
psi <- sqrt(h$counts / normFactor)  # normalized wave function # nolint

sum(psi ^ 2) * binWidth  # should be approximately 1 # nolint

psiAnalytical <- exp(-(h$mids ^ 2) / 2) # Obtaining the analytical wavefunction #nolint

psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) # Normalising the analytical wavefunction #nolint

sum(psiAnalytical ^ 2) * binWidth

ggplot(data.frame(x = h$mids, psi = psi), aes(x = x, y = psi)) +
  geom_point() +
  geom_line(aes(x = h$mids, y = psiAnalytical, color = "blue")) +
  labs(x = "Position", y = "Psi") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))

correlationDQHO <- read.csv("csv/correlationDirichletQHO.csv") #nolint

head(correlationDQHO)

L <- 5000 #nolint

t <- seq(0, L, length.out = L)

correlator <- function(t) {
  0.5 * exp(-t) + 0.5 * exp(t - (L - 1))
}

gAnalytic <- data.frame( #nolint
  t = t,
  G = sapply(0:(L - 1), correlator)
)

ggplot(correlationDQHO, aes(x = as.numeric(row.names(correlationDQHO)), y = Correlation)) + #nolint
  geom_point() +
  geom_line(data = gAnalytic, aes(x = t, y = G), color = "red") +
  labs(x = "Distance", y = "Correlation") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))

min(correlationDQHO$Correlation)
which(correlationDQHO$Correlation == min(correlationDQHO$Correlation))

## Double Well Potential ##

thermDDWP <- read.csv("csv/E0ThermalisationDirichletDWP.csv") #nolint

ggplot(thermDDWP, aes(x = as.numeric(row.names(thermDDWP)), y = E0)) +
  geom_point()


E0evolutionDDWP <- read.csv("csv/E0EvolutionDirichletDWP.csv") #nolint
ggplot(E0evolutionDDWP, aes(x = as.numeric(row.names(E0evolutionDDWP)), y = E0)) + #nolint
  geom_point() + geom_line()

waveFunctionDDWP <- read.csv("csv/waveFunctionDirichletDWP.csv") #nolint
x <- waveFunctionDDWP$Position

bins <- 50

hist(x, breaks = bins, main = "Histogram of positions",
     xlab = "position", ylab = "count")
h <- hist(x, breaks = 200, plot = FALSE)

measures <- sum(h$counts)  # total counts which equals measures * lattice size
positionRange <- range(h$mids)[2] - range(h$mids)[1]  # range of positions # nolint
binWidth <- positionRange / bins  # width of each bin # nolint
normFactor <- measures * binWidth  # normalization factor # nolint
psi <- sqrt(h$counts / normFactor)  # normalized wave function # nolint

sum(psi ^ 2) * binWidth  # should be approximately 1 # nolint

psiAnalytical <- exp(-(h$mids ^ 2) / 2) # Obtaining the analytical wavefunction #nolint

psiAnalytical <- psiAnalytical / sqrt(sum(psiAnalytical ^ 2) * binWidth) # Normalising the analytical wavefunction #nolint

sum(psiAnalytical ^ 2) * binWidth

ggplot(data.frame(x = h$mids, psi = psi), aes(x = x, y = psi)) +
  geom_point() +
  geom_line(aes(x = h$mids, y = psiAnalytical, color = "blue")) +
  labs(x = "Position", y = "Psi") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))

correlationDDWP <- read.csv("csv/correlationDirichletDWP.csv") #nolint

head(correlationDDWP)

L <- 5000 #nolint

t <- seq(0, L, length.out = L)

correlator <- function(t) {
  0.5 * exp(-t) + 0.5 * exp(t - (L - 1))
}

gAnalytic <- data.frame( #nolint
  t = t,
  G = sapply(0:(L - 1), correlator)
)

ggplot(correlationDDWP, aes(x = as.numeric(row.names(correlationDDWP)), y = Correlation)) + #nolint
  geom_point() +
  geom_line(data = gAnalytic, aes(x = t, y = G), color = "red") +
  labs(x = "Distance", y = "Correlation") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))

min(correlationDDWP$Correlation)
which(correlationDDWP$Correlation == min(correlationDDWP$Correlation))