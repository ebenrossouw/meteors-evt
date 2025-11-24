library(ReIns)
library(ismev)

options(scipen=999)
par.default = par(no.readonly = TRUE)

# GPD Functions
H = function(y, g, s) {
    1 - (1 + g * y / s)^(-g^-1)
}

Q = function(p, g, s) {
    s/g * ((1 - p) ^ -g - 1)
}

h = function(y, g, s) {
    (1 + g * y / s)^-(1 + 1 / g) / s
}

l = function(y, g, s) {
    # (4.10) in S. Coles
    k = length(y)
    -k * log(s) - (1 + 1/g) * sum(log(1 + g * y / s))
}


## Data

fireballs = read.csv("data/cneos/cneos_fireball_data.csv")
colnames(fireballs) = c("time", "lat", "long", "alt", "vx", "vy", "vz", "E0", "E")
fireball_maxima = read.csv("data/cneos/cneos_fireball_12m.csv")

n = nrow(fireballs)
x = sort(fireballs$E)
c = fireball_maxima$count[-33][-1]


## Fit GPD

gpd = GPDmle(x)
k = gpd$k
g = gpd$gamma
s = gpd$sigma
u = x[n - k]

ss = s - g * u
e = MeanExcess(x, plot=FALSE)$e

e.se = numeric(n-1)

for (i in k) {
    y = x[x > u[i]] - u[i]
    e.se[i] = sd(y) / sqrt(length(y))
}

e.ci = cbind(e - e.se * qnorm(0.975), e + e.se * qnorm(0.975))


## Threshold Selection

# (4.9) in S. Coles
a = s / (1 - g)
b = g / (1 - g)

{
    plot(
        u, e,
        type="l",
        lwd=1.5,
        ylim=c(min(e.ci, na.rm = TRUE), max(e.ci, na.rm = TRUE)),
        main="Mean Residual Life Plot",
        xlab="Threshold",
        ylab="Mean Excess"
    )
    abline(h=0)
    
    lines(u, e.ci[,1], lwd=1.5, col="gray60")
    lines(u, e.ci[,2], lwd=1.5, col="gray60")
    
    j = 90
    abline(a[j], b[j], lty=2)
}

{
    i = 1:200
    
    plot(
        k[i],
        b[i],
        type="l",
        lwd=1.5,
        main="Estimated Slope of Mean Residual Life Plot",
        xlab="Threshold",
        ylab="Slope of Mean Excess"
    )
    
    abline(h=0, v=0)
    abline(v=i[i %% 10 == 0], lty="dashed", col="lightgray")
}

{
    lim = c(0, 441)
    # lim = c(0, 60)

    plot(
        NA,
        ylim = lim,
        xlim = lim,
        xlab="Fitted Quantiles",
        ylab="Sample Quantiles"
    )
    title(main="QQ Plot for Fitted GPD on Impact Energy \nExceedances")
    mtext("Over Candidate Thresholds k=100, 90, 80, 70", side = 3, line = 0, adj = 0.5)
    abline(0, 1)
    
    # col = c("gray75", "gray50", "gray0")
    # i = c(100, 90, 70)
    col = c("gray80", "gray60", "gray40", "gray0")
    i = c(100, 90, 80, 70)
    
    for (j in 1:4) {
        y = x[x > u[i[j]]] - u[i[j]]
        N = length(y)
        q = Q(1:N / (N+1), g[i[j]], s[i[j]])
        
        points(q, y, pch=4, col=col[j], lwd=2)
    }

    legend(
        "bottomright",
        legend = paste("k =", k[i]),
        col=col,
        lty=0,
        pch=4,
        lwd=2
    )
}

## Inference

i = 90
y = x[x > u[i]] - u[i]
z = length(y) / n
# N = length(y)

# V = matrix(0, 3, 3)
# V[1:2, 1:2] = gpd.fit(x, u[i], siginit=s[i], shinit=g[i])$cov
# V[3,3] = z * (1 - z) / n

se = gpd.fit(x, u[i], siginit=s[i], shinit=g[i])$se

s.se = se[1]
g.se = se[2]
z.se = sqrt(z * (1 - z) / n)

round(rbind(
    c(s[i], s.se, s[i] + c(-1, 1) * s.se * qnorm(0.975)),
    c(g[i], g.se, g[i] + c(-1, 1) * g.se * qnorm(0.975)),
    c(z,    z.se, z    + c(-1, 1) * z.se * qnorm(0.975))
), 6)

q = as.integer((441 / 2):(441 * 3 + 1))
p = z * (1 - H(q - u[i], g[i], s[i]))
m = 1 / p
r = m / mean(c)

# d = cbind(
#     # d/ds p
#     z ^ ((2 + g[i]) / (1 + g[i])) * (q - u[i]) / s[i]^2 * p ^ (1 + g[i]),
#     # d/dg p
#     -(q - u[i]) / s[i] / g[i] * z ^ (1 / (1 + g[i])) * p ^ (1 + g[i])
#         + p * log(1 + g[i] * (q - u[i]) / s[i]) / g[i],
#     # d/dz p
#     p / z
# )
# 
# p.var = diag(d %*% V %*% t(d))
# p.se = sqrt(p.var)
# p.ci2 = cbind(p - p.se * qnorm(0.975), p + p.se * qnorm(0.975))

B = length(p)
p.ci = matrix(NA, B, 2)
p.l = l(y, g[i], s[i])

lp = function(y, p, q, z) {
    optimise(
        function(g) {
            s = g * q * ((p / z) ^ -g - 1)^-1
            l(y, g, s)
        },
        c(0.001, 2),
        maximum = TRUE
    )$objective
}

D = function(p, q) {
    p.pl = lp(y, p, q, z)
    D = 2 * (p.l - p.pl)
    D - qchisq(0.95, 1)
}

for (j in 1:B) {
    p.ci[j,1] = uniroot(D, c(1e-8, p[j]), q=q[j] - u[i], tol=1e-8)$root
    p.ci[j,2] = uniroot(D, c(p[j], 0.01), q=q[j] - u[i], tol=1e-8)$root
}

m.ci = 1 / p.ci[,2:1]
r.ci = m.ci / mean(c)


j = q %in% round(441 * 1:6 / 2)
round(cbind(
    q[j],
    p[j], p.ci[j,],
    m[j], m.ci[j,],
    r[j], r.ci[j,]
), 8)

{
    par(mar=c(4, 5, 3, 2) + 0.1)
    
    plot(
        q,
        p,
        type="l",
        ylim=c(min(p.ci), max(p.ci)),
        yaxt="n",
        ylab="",
        xlab="Impact Energy (kt)",
        lwd=2,
        main="Estimated GPD Exceedance Probabilities"
    )
    axis(2, las=2)
    title(ylab="Exceedance Probability", line=4)
    
    abline(h=-4:20 / 1e4, lty=2, col="lightgray")
    abline(h=0)
    
    lines(q, p.ci[,1], col="gray60", lwd=2)
    lines(q, p.ci[,2], col="gray60", lwd=2)
    
    j = which(q == 441)
    points(q[j], p[j], pch=4, cex=2)
    
    par(par.default)
}

{
    par(mar=c(4, 5, 3, 2) + 0.1)
    
    plot(
        q,
        r,
        type="l",
        ylim=c(min(r.ci), max(r.ci)),
        yaxt="n",
        ylab="",
        xlab="Impact Energy (kt)",
        lwd=1.5,
        main="Estimated GPD Return Periods in Years",
        log="y"
    )
    axis(2, las=2)
    title(ylab="Return Period (Years)", line=4)
    # abline(h=0:8 * 1e4 / 2, lty=2, col="lightgray")
    # abline(h=0:10 * 1e2, lty=2, col="lightgray")
    yticks=c(10 ^ (1:5), 5 * 10 ^ (1:5))
    abline(h=1*10 ^ (1:5), lty=2, col="gray80")
    abline(h=5*10 ^ (1:5), lty=2, col="gray80")
    abline(h=2*10 ^ (1:5), lty=2, col="gray80")
    
    lines(q, r.ci[,1], col="gray60", lwd=1.5)
    lines(q, r.ci[,2], col="gray60", lwd=1.5)
    
    j = which(q == 441)
    points(q[j], r[j], pch=4, cex=2)
    
    par(par.default)
}
























