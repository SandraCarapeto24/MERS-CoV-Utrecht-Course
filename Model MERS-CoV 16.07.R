#MERS-CoV2 Model 
#-------------
#Vaccination Combinations
vaccination_rates <- seq(0, 100, by = 20)

#Proportions of Vaccination possible combinations
combination_matrix <- outer(
  vaccination_rates,
  vaccination_rates,
  FUN = function(h, c) paste0(h, "%-", c, "%")
)

# Add row and column names for clarity
rownames(combination_matrix) <- paste0("Human_", vaccination_rates, "%")
colnames(combination_matrix) <- paste0("Camel_", vaccination_rates, "%")

# Print matrix
print(combination_matrix)
#------------------------

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

#Parameters 
pars <- c(betacc = 0.2, #beta = transmission rate
                 omegac = 0,
                 sigmac = 0.1429, # Latent period for camels 7 days
                 gammac = 0.0047, #gamma = recovery rate
          betahh = 0.2, #beta = transmission rate
          betahc = 0.2226,
          omegah = 0,
          sigmah = 0.1429, # Latent period for humans 7 days
          gammah = 0.0759, #gamma = recovery rate
          muh = 0.0019 # case fatality rate for humans 
)

#initial values for camels
N0c = 1800000 #population size for camels
pc <- 0
V0c = pc * N0c #initial fraction vaccinated for camels
Y0c = 0 #initial fraction infected for camels
Z0c = 0 #initial fraction recovered for camels
W0c = 0 #initial fraction latent for camels
X0c = (1-pc) * N0c - Y0c - W0c - Z0c

#initial values for humans
N0h = 33260000 #population size for humans
ph <- 0
V0h = ph * N0h #initial fraction vaccinated for humans
Y0h = 1 #initial fraction infected for humans
Z0h = 0 #initial fraction recovered for humans
W0h = 0 #initial fraction latent for humans
X0h = (1-ph) * N0h - Y0h - W0h - Z0h
init <- c(Xc = X0c,Vc = V0c,Wc = W0c, Yc = Y0c, Zc = Z0c, Xh = X0h,Vh = V0h,Wh = W0h, Yh = Y0h, Zh = Z0h)

#Data storage time
dt = 0.5#timestep for storing data
sim_duration = 1500 #length of the simulation
times <- seq(0, sim_duration, by = dt)

#Solve the ordinary differential equations
ode.out <- ode(init, times, sveir,pars)

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


#plot X against Y
plot(x = ode.out[,2], y = ode.out[,3], type = "l", xlab ="X", ylab ="Y")
s <- seq(length(ode.out[,2])-1)
arrows(x0 = ode.out[s,2], y0 = ode.out[s,3], x1 = ode.out[s + 1,2], y1 =
         ode.out[s + 1 ,3],length = .075)
