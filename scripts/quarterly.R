library(ReIns)
library(ismev)

fireballs = read.csv("data/cneos/cneos_fireball_quarterly.csv")

z = exp(fireballs$logE)

# Check Stationarity
acf(z)

# Fit GEV
gev = gev.fit(z)
cbind(gev$mle, gev$se)

gev.diag(gev)

G = function(z) {
    mu = gev$mle[1]
    sigma = gev$mle[2]
    gamma = gev$mle[3]
    
    exp(-(1 + gamma * (z - mu) / sigma) ^ (-1 / gamma))
}

1 - G(4:10 * 100)
1 - G(441)

1/(1 - G(4:10 * 100))
1/(1 - G(441))

q = 400:1000
plot(q, 1 - G(q), type="l")
plot(q, 1/(1 - G(q)), type="l")

