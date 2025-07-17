
library(deSolve)
########################################################################

#Ordinary differential Equations
sveir <- function(Time, State, pars){
  with(as.list(c(State,pars)),{
    Nc <- Xc + Vc + Wc + Yc + Zc  # Total camel population
    dXc <- - 1 * betacc * Xc * Yc/Nc + omegac * Vc 
    dVc <- -omegac * Vc 
    dWc <- betacc * Xc * Yc / Nc - sigmac * Wc
    dYc <- sigmac * Wc - gammac * Yc
    dZc <- gammac * Yc 
    Nh <- Xh + Vh + Wh + Yh + Zh  # Total camel population
    dXh <- - 1 * betahh * Xh * Yh/Nh - 1 * betahc * Xh * Yc/Nc + omegah * Vh 
    dVh <- -omegah * Vh 
    dWh <- betahh * Xh * Yh / Nh + betahc * Xh * Yc / Nc - sigmah * Wh
    dYh <- sigmah * Wh - gammah * Yh - muh * Yh
    dZh <- gammah * Yh 
    return(list(c(dXc, dVc, dWc, dYc, dZc, dXh, dVh, dWh, dYh, dZh)))
  })
}

param_table <- read.csv("Data.csv", sep = ";", stringsAsFactors = FALSE)

run_index <- 1

if (run_index < 1 || run_index > nrow(param_table)) {
  stop("Invalid run index. Please choose a number between 1 and ", nrow(param_table))
}

row <- param_table[run_index, ]

cat("\n--- Running scenario", run_index, "---\n")

# Extract and construct parameter vector for the model
pars <- c(
  betacc = row$betacc,
  omegac = row$omegac,
  sigmac = row$sigmac,
  gammac = row$gammac,
  betahh = row$betahh,
  betahc = row$betahc,
  omegah = row$omegah,
  sigmah = row$sigmah,
  gammah = row$gammah,
  muh = row$muh
)

# Initial values from selected row
N0c <- row$N0c
pc  <- row$pc
Y0c <- row$Y0c
Z0c <- row$Z0c
W0c <- row$W0c
V0c <- pc * N0c
X0c <- (1 - pc) * N0c - Y0c - Z0c - W0c

N0h <- row$N0h
ph  <- row$ph
Y0h <- row$Y0h
Z0h <- row$Z0h
W0h <- row$W0h
V0h <- ph * N0h
X0h <- (1 - ph) * N0h - Y0h - Z0h - W0h

init <- c(Xc = X0c, Vc = V0c, Wc = W0c, Yc = Y0c, Zc = Z0c,
          Xh = X0h, Vh = V0h, Wh = W0h, Yh = Y0h, Zh = Z0h)

# Print parameter values
cat("Parameters used in the model run:\n")
print(pars)

# Print initial conditions
cat("\nInitial values used:\n")
print(init)

# Time settings
dt <- 0.5
sim_duration <- 1500
times <- seq(0, sim_duration, by = dt)

# Solve ODEs
ode.out <- ode(init, times, sveir, pars)

# Set up a 2x1 layout with extra space for legend
layout(matrix(c(1, 2, 3, 4), nrow = 4), heights = c(4, 1, 4, 1))

# Increase margins and font
par(mar = c(5, 5, 4, 2), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.8)

## --- Camel Dynamics Plot ---
plot(ode.out[, 1], ode.out[, 2], type = "l", col = "blue", lwd = 3,
     ylim = c(0, N0c),
     xlab = "Time (days)", ylab = "Number of Camels",
     main = "Camel Dynamics")
lines(ode.out[, 1], ode.out[, 3], col = "purple", lwd = 3)
lines(ode.out[, 1], ode.out[, 4], col = "orange", lwd = 3)
lines(ode.out[, 1], ode.out[, 5], col = "red", lwd = 3)
lines(ode.out[, 1], ode.out[, 6], col = "green", lwd = 3)

## --- Camel Legend (in its own space) ---
par(mar = c(0, 0, 0, 0))  # remove margins for legend
plot.new()
legend("center", legend = c("Xc", "Vc", "Wc", "Yc", "Zc"),
       col = c("blue", "purple", "orange", "red", "green"),
       lwd = 3, horiz = TRUE, cex = 1.2, bty = "n", text.font = 2)

## --- Human Dynamics Plot ---
par(mar = c(5, 5, 4, 2))  # reset margins
plot(ode.out[, 1], ode.out[, 7] / 1e6, type = "l", col = "blue", lwd = 3,
     ylim = c(0, N0h / 1e6),
     xlab = "Time (days)", ylab = "Number of Humans (millions)",
     main = "Human Dynamics")
lines(ode.out[, 1], ode.out[, 8] / 1e6, col = "purple", lwd = 3)
lines(ode.out[, 1], ode.out[, 9] / 1e6, col = "orange", lwd = 3)
lines(ode.out[, 1], ode.out[, 10] / 1e6, col = "red", lwd = 3)
lines(ode.out[, 1], ode.out[, 11] / 1e6, col = "green", lwd = 3)

## --- Human Legend (own space) ---
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Xh", "Vh", "Wh", "Yh", "Zh"),
       col = c("blue", "purple", "orange", "red", "green"),
       lwd = 3, horiz = TRUE, cex = 1.2, bty = "n", text.font = 2)

