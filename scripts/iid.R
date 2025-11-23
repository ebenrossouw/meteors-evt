library(ReIns)
options(scipen=999)
par.default = par(no.readonly = TRUE)

# Data

fireballs = read.csv("data/cneos/cneos_fireball_data.csv")
colnames(fireballs) = c("time", "lat", "long", "alt", "vx", "vy", "vz", "E0", "E")
fireball_maxima = read.csv("data/cneos/cneos_fireball_12m.csv")

n = nrow(fireballs)
x = sort(fireballs$E)
c = fireball_maxima$count[-33][-1]

# QQ-Plots

ParetoQQ(x)
ParetoQQ_der(x, ylim=c(0, 21))  # Note y-axis scale

genQQ(x, gamma=hill$gamma)

{
    par(mfrow=c(4, 2), mar=c(4, 4, 3, 2) + 0.1, pch=4)
    
    ParetoQQ(x)
    abline(a=-h*qexp((n-k)/(n+1))+log(x[n-k]), b=h, lwd=1.5, col="blue")
    abline(v=qexp((n-k)/(n+1)), h=log(x[n-k]), lty=2, col="gray")
    ParetoQQ_der(x)
    abline(h=h, lwd=1.5, col="blue")
    
    LognormalQQ(x)
    LognormalQQ_der(x)
    
    ExpQQ(x)
    MeanExcess(x)
    
    WeibullQQ(x)
    WeibullQQ_der(x)
    
    par(par.default)
}


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

round(rbind(
    c(h,  abias, (avar),    ci),
    c(g, 0,      (avar.2o), ci.2o)
), 6)

{
    j = hill$k[1:n]
    plot(
        NA,
        ylim=c(0.5, 1.5),
        xlim=c(1, n),
        ylab="EVI",
        xlab="k",
        main="Pareto-Type Estimators of EVI"
    )
    
    # lines(j, epd$gamma[j], lwd=1.5, col="green")
    lines(j, hill$gamma[j], lwd=1.5, col="blue")
    lines(j, hill2$gamma[j], lwd=1.5, col="red")
    
    abline(v=80, h=hill2$gamma[80], lwd=1.5, lty=2, col="red")
    abline(v=k,  h=h,               lwd=1.5, lty=2, col="blue")
    legend(
        "bottomright",
        c("Hill", "Hill (2nd Order)", "EPD"),
        col=c("blue","red", "green"),
        lwd=1.5
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

{
    ParetoQQ(x)
    
    e = qexp((n - k) / (n + 1))
    abline(a=log(x[n-k]) - h * e, b=h, lwd=1.5, col="blue")
    abline(v=e, h=log(x[n-k]), lty=2, col="gray")
    
    ParetoQQ_der(x)
    abline(h=h, lwd=1.5, col="blue")
    abline(v=log(x[n-80]), lwd=1.5, col="red")
}

# Inference

q = as.integer((441 / 2):(441 * 3 + 1))
m = length(q)
p = numeric(m)
l = numeric(m)
u = numeric(m)

for (i in 1:m) {
    weissman = Weissman.p(x, hill$gamma, q[i])
    p[i] = weissman$P[k]
    
    d = z * (k * (1 + h ^ -2 * (log(q[i] / x[n-k])) ^ 2)) ^ -0.5
    l[i] = p[i] / (1 + d)
    u[i] = p[i] / (1 - d)
}

{
    par(mar=c(4, 5, 3, 2) + 0.1)

    plot(
        q,
        p,
        type="l",
        yaxt="n",
        ylab="",
        xlab="Impact Energy (kt)",
        lwd=1.0,
        main="Estimated Weissman Exceedance Probabilities"
    )
    axis(2, las=2)
    title(ylab="Exceedance Probability", line=4)
    
    lines(q, l, lty="dashed")
    lines(q, u, lty="dashed")
    
    i = which(q == 441)
    points(q[i], p[i], pch=3, cex=2)
    # abline(v=q[i], col="gray")
    # abline(v=q[i], h=p[i], lty=2, col="gray")
    abline(h=1:8 / 1e4, lty=2, col="lightgray")
    
    par(par.default)
}

{
    par(mar=c(4, 5, 3, 2) + 0.1)
    plot(
        q,
        1/p,
        type="l",
        xlab="Impact Energy (kt)",
        yaxt="n",
        ylab="",
        # ylab="Return Period",
        main="Estimated Weissman Return Period"
    )
    axis(2, las=2)
    title(ylab="Return Period", line=4)
    lines(q, 1/u, lty="dashed")
    lines(q, 1/l, lty="dashed")
    
    # j = which(q == 441)
    # abline(v=q[j], h=1/p[j], lty=2, col="gray")
    # par(par.default)
    i = which(q == 441)
    points(q[i], 1/p[i], pch=3, cex=2)
    # abline(v=q[i], col="gray")
    # abline(v=q[i], h=p[i], lty=2, col="gray")
    abline(h=1:8 * 1000, lty=2, col="lightgray")
}

# i = c(0:2 * 441 + 42, 0:10 * 100 + 1)
# i = round(220.5 * 2:6) - 400 + 1
i = q %in% round(441 * 1:6 / 2)
round(cbind(q[i], p[i], l[i], u[i], 1/p[i], 1/u[i], 1/l[i]), 6)

# ppareto(q[i], 1/g, lower.tail=FALSE) # TODO


## 2nd Order

k = 80
g = hill2$gamma
b = hill2$b
beta = hill2$beta

v = g[k] ^ 2 * ((1 + beta[k]) / beta[k]) ^ 4

q = as.integer((441 / 2):(441 * 3 + 1))
B = length(q)
p = numeric(B)

for (j in 1:B) {
    f = function(p, q) Quant.2oQV(x, g, b, beta, p)$Q[k] - q
    p[j] = uniroot(f, c(1e-8, 0.1), q=q[j], tol=1e-8)$root
}

m = 1 / p
r = m / mean(c)

p.ci = matrix(NA, B, 2)

for (j in 1:B) {
    # d = z * (k * (1 + (((1 + beta[k]) / beta[k]) ^ 2 * g[k]) ^ -2 * (log(q[j] / x[n-k])) ^ 2)) ^ -0.5
    d = z * (k * (1 + (1 / v) * (log(q[j] / x[n-k])) ^ 2)) ^ -0.5
    p.ci[j,1] = p[j] / (1 + d)
    p.ci[j,2] = p[j] / (1 - d)
}

m.ci = 1 / p.ci
r.ci = m.ci / mean(c)


{
    par(mar=c(4, 5, 3, 2) + 0.1)
    
    plot(
        NA,
        ylim=c(min(p.ci), max(p.ci)),
        ylab="",
        yaxt="n",
        xlim=c(min(q), max(q)),
        xlab="Impact Energy (kt)",
        main="Estimated Exceedance Probabilities"
    )
    axis(2, las=2)
    title(ylab="Exceedance Probability", line=4)
    
    abline(h=1:8 / 1e4 / 2, lty=2, col="lightgray")
    
    lines(q, p, lwd=1.5)
    lines(q, p.ci[,1], lwd=1.5, col="gray60")
    lines(q, p.ci[,2], lwd=1.5, col="gray60")
    
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






















