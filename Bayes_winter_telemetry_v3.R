
# library(IPMbook)
library(jagsUI)

setwd("C:/Users/David.Pavlacky/Documents/Monitoring/IPM/Winter_telemetry/Bayesian_estimation/Winter_survival")

# Load RData from Q. Latif

load("C:/Users/David.Pavlacky/Documents/Monitoring/IPM/Winter_telemetry/Bayesian_estimation/Winter_survival/Data_compiled_MissingCovsImputed.RData")

# ls()

# Grab encounter history data

data_bais <- as.data.frame(data.BAIS)

# write.csv(data_bais, "data_bais.csv")

# Extract site, year, encounter summary stats, and encounter histories

data_bais_2 <- data_bais[, c(1, 2, 4:6, 203:326)]
  names(data_bais_2) [1:5] <- c("Site", "Year", "f", "j", "k")

# Date when the nest is first detected (fi)
# Last day when nest i was encountered alive (ji)
# Day when nest i was visited for the last time (ki)
# Whether nest i was alive during this last visit (xi = 0 for a failed nest and xi = 1 for a successful nest)

# Specify individual fates
    
data_bais_2$x <- ifelse(apply(data_bais_2[, 6:129], 1, max, na.rm = T) > 1, 0, 1)
    
# Format factors

data_bais_2$Site <- as.factor(data_bais_2$Site)
data_bais_2$Year <- as.factor(data_bais_2$Year)
data_bais_2$Site_year <- as.factor(paste(data_bais_2$Site, data_bais_2$Year, sep = "-"))

# Specify indices for random effects

data_bais_2$year <- match(data_bais_2$Year, unique(data_bais_2$Year))
data_bais_2$site <- match(data_bais_2$Site, unique(data_bais_2$Site))
data_bais_2$site_yr <- match(data_bais_2$Site_year, unique(data_bais_2$Site_year))

# Extract summary stats, factors and indices 

bais <- data_bais_2[, c("f", "j", "k", "x", "Site", "Year", "Site_year", "year", "site", "site_yr")]

# Set k = j for surviving individuals never encountered again, right censor

bais$k <- ifelse(bais$x == 1, bais$j, bais$k)

# Manual correct data

bais[31, "j"] <- 104
bais[295, "j"] <- 63
bais[319, "j"] <- 78
bais[387, "j"] <- 52

# Remove empty time occasion 1

bais$f <- bais$f - 1
bais$j <- bais$j - 1
bais$k <- bais$k - 1

# Specify trap effect covariate for lower survival from 1 to 7 days after capture

bais$trap <- ifelse(bais$k - bais$f < 8, 0, 1)

# Implement nest survival model (Schaub and Kéry 2022)
# Schaub, M., and M. Kéry. 2022. Integrated population models: theory and ecological applications with R and JAGS. Academic Press, London, United Kingdom.
# Gives identical results to the S(.) nest survival model in RMark

# Identify failed broods

fail <- which(bais$x == 0)

# Create encounter histories

y <- matrix(NA, nrow = length(bais$f), ncol = max(bais$k))

for (i in 1:length(bais$f)) {
  y[i, bais$f[i]] <- 1
  y[i, bais$j[i]] <- 1
}

for (i in 1:length(fail)) {
  y[fail[i], bais$k[fail[i]]] <- 0
}

# Bundle data
jags.data <- with(bais, list(y = y, f = f, k = k, n.nest = nrow(y), n.day = max(k), n.site = length(unique(bais$Site)), n.year = length(unique(bais$Year)), n.site.yr = length(unique(bais$Site_year)), year = year, site = site, trap = trap, site_yr = site_yr, T = 107))
str(jags.data)

# Write JAGS model file

cat(file = "winter_survival.txt", "
    model {
    # Priors and linear models
    for (i in 1:n.nest) {
      for (t in f[i]:(k[i] - 1)) {
        # logit(phi[i, t]) <- alpha # S(.) identical to RMark
        # logit(phi[i, t]) <- alpha + beta * trap[i] + eta[year[i]] # Year random effect
        logit(phi[i, t]) <- alpha + beta * trap[i] + eta[site_yr[i]] # Site x year random effect
      } #t
    } #i

    # for (i in 1:n.year) { # Year random effect predictions  
    for (i in 1:n.site.yr) { # Site x year random effect predictions
      for (a in 1:n.nest) { 
        logit(phid[i, a]) <- alpha + beta + eta[i]
      }
    }
    
    alpha ~ dnorm(0, 0.001)
    beta ~ dnorm(0, 0.001)
    sigma_eta ~ dunif(0, 7)
    
    # for (t in 1:n.year) { # Year random effect
    for (t in 1:n.site.yr) { # Site x year random effect
      eta[t] ~ dnorm(0, pow(sigma_eta, -2))
    }
    
    # Likelihood
    for (i in 1:n.nest) {
      for (t in (f[i] + 1):k[i]) {
        y[i, t] ~ dbern(phi[i, t - 1] * y[i, t - 1])
      } #t
    } #i
    
    # Derived parameter: nest success, mid-season 90 day-survival
    # for (i in 1:n.year) { # Year random effect
    for (i in 1:n.site.yr) { # Site x year random effect
      nu[i] <- prod(phid[i, 18:T])
    }
    
    }
    ")

# Initial values
inits <- function(){list(alpha = runif(1, 4, 5), beta = runif(1, 0, 0.1))}

# Parameters monitored
parameters <- c("phi", "phid", "nu", "alpha", "beta", "eta", "sigma_eta")

# MCMC settings
ni <- 2000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS from R (ART 1 min) and check convergence
out <- jagsUI(jags.data, inits, parameters, "winter_survival.txt", n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = T)
print(out, 5)
# out$summary
# str(out)

out$Rhat$alpha
out$Rhat$beta
out$Rhat$eta
out$Rhat$sigma_eta

# dim(out$sims.list$phi)
# out$sims.list$phi[1, 4, ]

out$mean$alpha
out$sd$alpha

out$mean$beta
out$sd$beta

out$mean$eta
out$mean$sigma_eta

out$mean$nu
out$sd$nu
