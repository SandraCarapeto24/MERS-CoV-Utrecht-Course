
library(deSolve)

#Ordinary differential Equations
sveir_camels <- function(Time, State, pars_camels){
  with(as.list(c(State,pars_camels)),{
    Nc <- Xc + Vc + Wc + Yc + Zc  # Total camel population
    dXc <- - 1 * betacc * Xc * Yc/Nc + omegac * Vc 
    dVc <- -omegac * Vc 
    dWc <- betacc * Xc * Yc / Nc - sigmac * Wc
    dYc <- sigmac * Wc - gammac * Yc
    dZc <- gammac * Yc 
    return(list(c(dXc, dVc, dWc, dYc, dZc)))
  })
}
#Parameters for camels
pars_camels <- c(betacc = 0.0048, #beta = transmission rate
                 omegac = 0,
                 sigmac = 0.1429, # Latent period for camels 7 days
                 gammac = 0.0047 #gamma = recovery rate
)

#initial values for camels
N0c = 1800000 #population size for camels
pc <- 0.4
V0c = pc * N0c #initial fraction vaccinated for camels
Y0c = 1 #initial fraction infected for camels
Z0c = 0 #initial fraction recovered for camels
W0c = 0 #initial fraction latent for camels
X0c = (1-pc) * N0c - Y0c - W0c - Z0c
initc <- c(Xc = X0c,Vc = V0c,Wc = W0c, Yc = Y0c, Zc = Z0c) 


sveir_humans <- function(Time, State, pars_humans){
  with(as.list(c(State,pars_humans)),{
    Nh <- Xh + Vh + Wh + Yh + Zh  # Total camel population
    dXh <- - 1 * betahh * Xh * Yh/Nh + omegah * Vh 
    dVh <- -omegah * Vh 
    dWh <- betahh * Xh * Yh / Nh + betahc * Xh * Yc / Nc - sigmah * Wh
    dYh <- sigmah * Wh - gammah * Yh - muh * Yh
    dZh <- gammah * Yh 
    return(list(c(dXh, dVh, dWh, dYh, dZh)))
  })
}
#Parameters for humans
pars_humans <- c(betahh = 0.0896, #beta = transmission rate
                 betahc = 0.2226,
                 omegah = 0,
                 sigmah = 0.1429, # Latent period for humans 7 days
                 gammah = 0.0759, #gamma = recovery rate
                 muh = 0.0019 # case fatality rate for humans 
)

#initial values for humans
N0h = 33260000 #population size for humans
ph <- 0
V0h = ph * N0h #initial fraction vaccinated for humans
Y0h = 1 #initial fraction infected for humans
Z0h = 0 #initial fraction recovered for humans
W0h = 0 #initial fraction latent for humans
X0h = (1-ph) * N0h - Y0h - W0h - Z0h
inith <- c(Xh = X0h,Vh = V0h,Wh = W0h, Yh = Y0h, Zh = Z0h) 

#Solve the ordinary differential equations
ode.out_camels <- ode(initc, times, sveir_camels,pars_camels)

#Solve the ordinary differential equations
ode.out_humans <- ode(inith, times, sveir_humans,pars_humans)

#plot X, Y, and Z against time camels
plot(x = ode.out_camels[,1], y = ode.out_camels[,2], type = "l", xlab = "Time", ylab = "Count",
     col = "purple", lwd = 2, ylim = c(0, N0c))

lines(ode.out_camels[,1], ode.out_camels[,3], col = "red", lwd = 2)
lines(ode.out_camels[,1], ode.out_camels[,4], col = "palegreen1", lwd = 2)
lines(ode.out_camels[,1], ode.out_camels[,5], col = "steelblue", lwd = 2)
lines(ode.out_camels[,1], ode.out_camels[,6], col = "orange", lwd = 2)

legend("topright", legend = c("X", "V", "W", "Y", "Z"), 
       col = c("purple", "red", "palegreen1", "steelblue", "orange"),
       lwd = 2)

#plot X, Y, and Z against time humans
plot(x = ode.out_humans[,1], y = ode.out_humans[,2], type = "l", xlab = "Time", ylab = "Count",
     col = "purple", lwd = 2, ylim = c(0, N0c))

lines(ode.out_humans[,1], ode.out_humans[,3], col = "red", lwd = 2)
lines(ode.out_humans[,1], ode.out_humans[,4], col = "palegreen1", lwd = 2)
lines(ode.out_humans[,1], ode.out_humans[,5], col = "steelblue", lwd = 2)
lines(ode.out_humans[,1], ode.out_humans[,6], col = "orange", lwd = 2)

legend("topright", legend = c("X", "V", "W", "Y", "Z"), 
       col = c("purple", "red", "palegreen1", "steelblue", "orange"),
       lwd = 2)

