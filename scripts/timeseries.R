library(ReIns)
library(ismev)

par.default = par(no.readonly = TRUE)

# fb = read.csv("data/cneos/cneos_fireball_data.csv")
# 
# x = fb$Calculated.Total.Impact.Energy..kt.
# t = fb$Peak.Brightness.Date.Time..UT.
# # t = as.POSIXlt(t, format = "%Y-%m-%d %H:%M:%S")
# t = as.Date(t)
# t = t[t >= "1993-07-01"]


# cut.Date(t, "6 month")

# b = as.Date(fireballs$time)
# k = length(b)
# f = stepfun(b, c(NA, b))
# 
# r = mean(table(f(t)))
# 
# theta = function(u) {
#     mean(z > u) / (r * mean(x > u)) 
# }
# 
# v = quantile(x, 1:100/101)
# th = numeric(length(v))
# 
# for (i in 1:length(v)) {
#     th[i] = theta(v[i])
# }

# Q = function(p, gev) {
#     m = gev$mle[1]
#     s = gev$mle[2]
#     g = gev$mle[3]
#     
#     m - s/g * (1 - (-log(p))^-g)
# }
# 
# 
# for (m in c(1, 3, 6, 12, 18, 24)) {
#     path = sprintf("data/cneos/cneos_fireball_%dm.csv", m)
# 
#     fb = read.csv(path)
#     n = nrow(fb) - 1
#     mdl = gev.fit(fb$E[1:n])
#     
# 
#     print(path)
#     print(-mdl$nllh / n)
#     print(cor(Q(1:n / (n+1), mdl), sort(fb$E[1:n])))
# }

G = function(z, theta) {
    m = theta[1]
    s = theta[2]
    g = theta[3]
    exp(-(1 + g * (z - m) / s) ^ (-1 / g))
}

Q = function(p, theta) {
    m = theta[1]
    s = theta[2]
    g = theta[3]
    m - s/g * (1 - (-log(p))^-g)
}


# fireballs = read.csv("data/cneos/cneos_fireball_quarterly.csv")
fireballs = read.csv("data/cneos/cneos_fireball_6m.csv")

# drop first 3 and last 1
# z = fireballs$E[4:129]
z = fireballs$E[1:64]
n = length(z)

# Check Dependence
acf(z, main="Autocorrelation of Biannual Maximum Fireball Impact Energies")

# Fit GEV
gev = gev.fit(z)
# gev.diag(gev)


ci.l = gev$mle - qnorm(0.975) * gev$se;
ci.u = gev$mle + qnorm(0.975) * gev$se;

cbind(gev$mle, gev$se, ci.l, ci.u)
cor(Q(1:n / (n+1), gev$mle), sort(z))

{
    plot(
        Q(1:n / (n+1), gev$mle),
        sort(z),
        xlim=c(0, max(z)),
        # xlim=c(0, 59),
        # ylim=c(0, 59),
        main="QQ-Plot of GEV Fit on Biannual Fireballs",
        xlab="Theoretical Quantiles",
        ylab="Sample Quantiles"
    )
    abline(a = 0, b = 1)
}

{
    m = gev$mle[1]
    s = gev$mle[2]
    g = gev$mle[3]

    p = 1:1000/1001
    y = -log(1 - p)
    r = m - (s/g) * (1 - y^-g)
    
    plot(log(y), r)
}

q = as.integer((441 / 2):(441 * 3 + 1))
i = q %in% round(441 * 1:6 / 2)

1 - G(q[i], gev$mle)
1/(1 - G(q[i], gev$mle))

rbind(
    q[i],
    1 - G(q[i], gev$mle),
    1/(1 - G(q[i], gev$mle))
)

# Use T-squared to find simultaneous CIs for mean parameters */
# https://blogs.sas.com/content/iml/2016/12/07/simultaneous-confidence-intervals-mean.html
# k = 3
# f = qf(0.95, k, n-k)           # critical value of F(k, n-k) */
# T2 = k * (n - 1) / (n - k) * f  # Hotelling's T-squared is scaled F */
# ho.l = gev$mle - sqrt(T2) * gev$se;
# ho.u = gev$mle + sqrt(T2) * gev$se;

{
    par(mar=c(4, 5, 3, 2) + 0.1)
    
    plot(
        q,
        1 - G(q, gev$mle),
        type="l",
        yaxt="n",
        ylab="",
        # ylim=c(0, 0.005),
        xlab="Impact Energy (kt)",
        lwd=1.5,
        main="Estimated GEV Survival Function of Impact Energy Maxima per Half-Year"
    )
    axis(2, las=2)
    title(ylab="Exceedance Probability", line=4)
    
    # abline(v=441, h=1 - G(441), lty=2, col="gray")
    points(441, 1 - G(441, gev$mle), pch=3, cex=2)
    abline(h=0)
    abline(h=1:9 / 1000, lty=2, col="lightgray")
    
    par(par.default)
}

{
    par(mar=c(4, 5, 3, 2) + 0.1)
    
    plot(
        q,
        1/(1 - G(q, gev$mle)),
        type="l",
        yaxt="n",
        ylab="",
        # ylim=c(0, 0.005),
        xlab="Impact Energy (kt)",
        lwd=1.5,
        main="Estimated Return Periods of Biannual Maximum Fireball Impact Energy"
    )
    axis(2, las=2)
    title(ylab="Return Period", line=4)
    
    # abline(v=441, h=1 - G(441), lty=2, col="gray")
    points(441, 1/(1 - G(441, gev$mle)), pch=3, cex=2)
    abline(h=0)
    abline(h=1:10 * 100, lty=2, col="lightgray")
    
    par(par.default)
}
