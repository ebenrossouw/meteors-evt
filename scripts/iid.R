library(ReIns)

par.default = par(no.readonly = TRUE)

# Data
fireballs = read.csv("data/cneos/cneos_fireball_data.csv")
colnames(fireballs) = c("time", "lat", "long", "alt", "vx", "vy", "vz", "E0", "E")

n = nrow(fireballs)
x = sort(fireballs$E)
# x = log(x) + 4

# QQ-Plots
ParetoQQ(x)
ParetoQQ_der(x, ylim=c(0, 21))  # Note y-axis scale

LognormalQQ(x)
LognormalQQ_der(x, ylim=c(0, 21))

ExpQQ(x)
MeanExcess(x)

WeibullQQ(x)
WeibullQQ_der(x, ylim=c(0, 21))


# Estimators
hill = Hill(x)
hill2 = Hill.2oQV(x)
epd = EPD(x)

amse = hill2$gamma ^ 2 / hill2$k + (hill2$b / (1 + hill2$beta)) ^ 2
k = which.min(amse)

h = hill$gamma[k]
g = hill2$gamma[k]
b = hill2$b[k]
beta = hill2$beta[k]


# abias = b[k] * sum((1:k / (k+1)) ^ beta[k]) / k
abias = b / (1 + beta)
avar = h^2 / k
avar.2o = g^2 / k * (1 + 1/beta)^4

z = qnorm(0.975)
ci = h * c(1 + z / sqrt(k), 1 - z / sqrt(k))^-1
ci.2o = g * c(1 + z / sqrt(k), 1 - z / sqrt(k))^-1 * (1 + 1/beta)^-2

rbind(
    c(h,  abias, avar,    ci),
    c(g, 0,     avar.2o, ci.2o)
)

genQQ(x, gamma=hill$gamma)

{
    j = hill$k[1:n]
    plot(NA, ylim=c(0.5, 1.5), xlim=c(1, n))
    
    lines(j, hill$gamma[j], lwd=1.5, col="blue")
    lines(j, epd$gamma[j], lwd=1.5, col="red")
    lines(j, hill2$gamma[j], lwd=1.5, col="green")
    
    abline(v=k, h=h, lty=2)
    # legend("bottomright", c("Hill","EPD"), col=c("blue","orange"), lwd=1.5)
}

ghill = genHill(x, gamma=hill$gamma)
mme = Moment(x)
mle = GPDmle(x)

{
    j = ghill$k
    plot(NA, ylim=c(0.5, 1.5), xlim=c(0, n))
    
    lines(j, ghill$gamma[j], col="blue", lwd=1.5)
    lines(j, mme$gamma[j], col="red", lwd=1.5)
    lines(j, mle$gamma[j], col="green", lwd=1.5)
    
    abline(v=k, h=h, lty=2, col="gray")
}

ParetoQQ(x)
abline(a=-h*qexp((n-k)/(n+1))+log(x[n-k]), b=h, lwd=1.5, col="blue")
abline(v=qexp((n-k)/(n+1)), h=log(x[n-k]), lty=2, col="gray")
ParetoQQ_der(x)
abline(h=h, lwd=1.5, col="blue")


# Small Exceedences

q = 400:1000
m = length(q)
p = numeric(m)

for (i in 1:m) {
    weissman = Prob(x, hill$gamma, q[i])
    p[i] = weissman$P[k]
}

{
    plot(q, p, type="l", yaxt="n")
    axis(2, las=2)
    
    j = which(q == 441)
    abline(v=q[j], h=p[j], lty=2, col="gray")
}

i = c(400, 441, 500, 600, 700, 800, 900, 1000) - 400 + 1
cbind(q[i], p[i])


# ppareto(q[i], 1/g, lower.tail=FALSE) # TODO


# Return Periods




























