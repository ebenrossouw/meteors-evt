library(ReIns)
options(scipen=999)
par.default = par(no.readonly = TRUE)

## Data

fireballs = read.csv("data/cneos/cneos_fireball_data.csv")
colnames(fireballs) = c("time", "lat", "long", "alt", "vx", "vy", "vz", "E0", "E")
fireball_maxima = read.csv("data/cneos/cneos_fireball_12m.csv")

n = nrow(fireballs)
x = sort(fireballs$E)
c = fireball_maxima$count[-33][-1]


## Estimators

hill = Hill(x)
hill2 = Hill.2oQV(x)
epd = EPD(x)

amse = hill$gamma ^ 2 / hill2$k + (hill2$b / (1 + hill2$beta)) ^ 2
avar.2o = hill2$gamma ^ 2 / k * (1 + 1 / hill2$beta) ^ 4

{
    plot(NA, ylim=c(0, 0.15), xlim=c(0, n))
    lines(hill$k, amse, col="blue", lwd=1.5)
    lines(hill2$k, avar.2o, col="red", lwd=1.5)
}

min(amse, na.rm = TRUE)
which.min(amse)
min(avar.2o, na.rm = TRUE)
which.min(avar.2o)


{
    j = hill$k[1:n]
    plot(
        NA,
        ylim=c(0.5, 1.5),
        xlim=c(1, 400),
        ylab="EVI",
        xlab="k",
        main="Pareto-Type Estimators of the EVI"
    )
    
    lines(j, epd$gamma[j], lwd=2, col="gray80")
    lines(j, hill$gamma[j], lwd=2, col="blue")
    lines(j, hill2$gamma[j], lwd=2, col="red")
    
    segments(80,  hill2$gamma[80], -100, hill2$gamma[80], lwd=1.5, lty=2, col="red")
    segments(250, hill$gamma[250], -100, hill$gamma[250], lwd=1.5, lty=2, col="blue")
    segments(80,  0, 80, hill2$gamma[80],  lty=2, lwd=1.5, col="red")
    segments(250, 0, 250, hill$gamma[250], lty=2, lwd=1.5,  col="blue")
    
    legend(
        "bottomright",
        c("Hill", "Hill (2nd Order)", "EPD"),
        col=c("blue","red", "gray80"),
        lwd=2
    )
}

ghill = genHill(x, gamma=hill$gamma)
mme = Moment(x)
mle = GPDmle(x)

{
    j = ghill$k
    plot(
        NA,
        ylim=c(0.5, 1.5),
        xlim=c(0, n),
        ylab="EVI",
        xlab="k",
        main="Generalised Estimators of EVI"
    )
    
    lines(j, ghill$gamma[j], col="blue", lwd=1.5)
    lines(j, mle$gamma[j], col="green", lwd=1.5)
    lines(j, mme$gamma[j], col="red", lwd=1.5)
    
    abline(v=k, h=h, lty=2, col="gray")
    
    legend(
        "bottomright",
        c("Generalised Hill", "MLE", "MoM"),
        col=c("blue","green", "red"),
        lwd=1.5
    )
}


## QQ-Plots

ParetoQQ(x)
ParetoQQ_der(x, ylim=c(0, 21))  # Note y-axis scale

genQQ(x, gamma=hill$gamma, pch=16)

{
    par(mfrow=c(4, 2), mar=c(4, 4, 3, 2) + 0.1, pch=16)
    
    ParetoQQ(x)
    # abline(a=-h*qexp((n-k)/(n+1))+log(x[n-k]), b=h, lwd=1.5, col="blue")
    # abline(v=qexp((n-k)/(n+1)), h=log(x[n-k]), lty=2, col="gray")
    ParetoQQ_der(x)
    # abline(h=h, lwd=1.5, col="blue")
    
    LognormalQQ(x)
    LognormalQQ_der(x)
    
    ExpQQ(x)
    MeanExcess(x)
    
    WeibullQQ(x)
    WeibullQQ_der(x)
    
    par(par.default)
}

{
    par(mfrow=c(1, 2), mar=c(4, 4, 3, 2) + 0.1)

    k = 250
    h = hill$gamma[k]
    E = qexp((n - k) / (n + 1))
    
    ParetoQQ(x, xlim=c(E, 7), ylim=c(log(x[n-k]), log(max(x))))
    
    abline(a=log(x[n-k]) - h * E, b=h, lwd=2, col="blue")
    abline(v=E, h=log(x[n-k]), lty=2, col="gray")
    
    ParetoQQ_der(x)
    abline(h=h, lwd=1.5, col="blue")
    # abline(v=log(x[n-80]), lwd=1.5, col="red")
    
    par(par.default)
}


## Inference (1st Order)

k = 250
h = hill$gamma[k]

h.abias = hill2$b[k] / (1 + hill2$beta[k])
h.avar = h ^ 2 / k
z = qnorm(0.975)
h.ci = h * c(1 + z / sqrt(k), 1 - z / sqrt(k))^-1

round(c(h,  h.abias, h.avar, h.ci), 6)

q = as.integer((441 / 2):(441 * 3 + 1))
B = length(q)
p = numeric(B)
p.ci = matrix(NA, B, 2)

for (i in 1:B) {
    p[i] = Weissman.p(x, hill$gamma, q[i])$P[k]
    
    d = z * (k * (1 + h ^ -2 * (log(q[i] / x[n-k])) ^ 2)) ^ -0.5
    p.ci[i,1] = p[i] / (1 + d)
    p.ci[i,2] = p[i] / (1 - d)
}

m = 1 / p
r = m / mean(c)
m.ci = 1 / p.ci[,2:1]
r.ci = m.ci / mean(c)

i = q %in% round(441 * 1:6 / 2)
round(cbind(
    q[i],
    p[i], p.ci[i,],
    m[i], m.ci[i,],
    r[i], r.ci[i,]
), 6)

{
    par(mar=c(4, 5, 3, 2) + 0.1)

    plot(
        q,
        p,
        type="l",
        yaxt="n",
        ylab="",
        xlab="Impact Energy (kt)",
        lwd=1.5,
        main="Estimated Weissman Exceedance Probabilities"
    )
    axis(2, las=2)
    title(ylab="Exceedance Probability", line=4)
    abline(h=1:8 / 1e4, lty=2, col="lightgray")
    
    lines(q, p.ci[,1], lwd=1.5, col="gray60")
    lines(q, p.ci[,2], lwd=1.5, col="gray60")
    
    i = which(q == 441)
    points(q[i], p[i], pch=3, cex=2)
    
    par(par.default)
}

{
    par(mar=c(4, 5, 3, 2) + 0.1)
    plot(
        q,
        r,
        type="l",
        xlab="Impact Energy (kt)",
        yaxt="n",
        ylab="",
        # ylab="Return Period",
        main="Estimated Weissman Return Period in Years"
    )
    axis(2, las=2)
    title(ylab="Return Period (Years)", line=4)
    abline(h=1:8 * 25, lty=2, col="lightgray")
    
    lines(q, r.ci[,1], lwd=1.5, col="gray60")
    lines(q, r.ci[,2], lwd=1.5, col="gray60")
    
    i = which(q == 441)
    points(q[i], r[i], pch=3, cex=2)
}


## Inference (2nd Order)

k = 80
g = hill2$gamma
b = hill2$b
beta = hill2$beta

# v = g[k] ^ 2 * ((1 + beta[k]) / beta[k]) ^ 4
g.avar = g ^ 2 / k * ((1 + beta) / beta) ^ 4
z = qnorm(0.975)
g.ci = g * (1 + 1 / beta) ^ -2 * c(1 + z / sqrt(k), 1 - z / sqrt(k)) ^ -1

round(c(g[k], g.avar[k], g.ci[k]), 6)
b[k]
beta[k]


q = as.integer((441 / 2):(441 * 3 + 1))
B = length(q)
p = numeric(B)
p.ci = matrix(NA, B, 2)

for (j in 1:B) {
    f = function(p, q) Quant.2oQV(x, g, b, beta, p)$Q[k] - q
    # p[j] = uniroot(f, c(1e-8, 0.1), q=q[j], tol=1e-8)$root

    d = z * (k * (1 + ((1 + beta[k]) / beta[k]) ^ -4 * g[k] ^ -2 * (log(q[j] / x[n-k])) ^ 2)) ^ -0.5
    p.ci[j,1] = p[j] / (1 + d)
    p.ci[j,2] = p[j] / (1 - d)
}

m = 1 / p
r = m / mean(c)
m.ci = 1 / p.ci[,2:1]
r.ci = m.ci / mean(c)

i = q %in% round(441 * 1:6 / 2)
round(cbind(
    q[i],
    p[i], p.ci[i,],
    m[i], m.ci[i,],
    r[i], r.ci[i,]
), 6)


{
    par(mar=c(4, 5, 3, 2) + 0.1)
    
    plot(
        NA,
        ylim=c(min(p.ci), max(p.ci)),
        ylab="",
        yaxt="n",
        xlim=c(min(q), max(q)),
        xlab="Impact Energy (kt)",
        main="Estimated 2nd Order Exceedance Probabilities"
    )
    axis(2, las=2)
    title(ylab="Exceedance Probability", line=4)
    
    abline(h=1:8 / 1e4 / 2, lty=2, col="lightgray")
    
    lines(q, p.ci[,1], lwd=1.5, col="gray60")
    lines(q, p.ci[,2], lwd=1.5, col="gray60")
    lines(q, p, lwd=1.5)
    
    j = which(q == 441)
    points(q[j], p[j], pch=3, cex=2)
    
    par(par.default)
}

{
    par(mar=c(4, 5, 3, 2) + 0.1)
    
    plot(
        NA,
        ylim=c(min(r.ci), max(r.ci)),
        ylab="",
        yaxt="n",
        xlim=c(min(q), max(q)),
        xlab="Impact Energy (kt)",
        main="Estimated Return Periods"
    )
    axis(2, las=2)
    title(ylab="Return Period", line=4)
    
    abline(h=1:16 * 100, lty=2, col="lightgray")
    
    lines(q, r, lwd=1.5)
    lines(q, r.ci[,1], lwd=1.5, col="gray60")
    lines(q, r.ci[,2], lwd=1.5, col="gray60")
    
    j = which(q == 441)
    points(q[j], r[j], pch=3, cex=2)
    
    par(par.default)
}






















