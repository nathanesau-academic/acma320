# Acma 320 - Spring 2014
# Assignment 3 
# Nathan Esau

a <- 0.5
x <- 10
n <- 5
omega <- 20
radix <- 1e+05

tpx <- function(t, x, a, omega) {
  (1 - (x+t)/omega)^a / (1 - x/omega)^a
}

tqx <- function(t, x, a, omega) {
  1 - tpx(t, x, a, omega)
}

createLifeTable <- function(x, radix, a, omega) {
  
  k = 0:(omega - x - 1)
  lx = sapply(k, function(t) tpx(t, x, a, omega)) * radix
  dx = c(-diff(lx), tail(lx, 1))
  qx <- sapply(k, function(t) tqx(1, x + t, a, omega))
  px <- 1 - qx
  
  data.frame(x = x + k, lx = lx, dx = dx, qx = qx, px = px)
}

lifeTable <- createLifeTable(x, radix, a, omega)

plot(function(t) tpx(t, x, a, omega), 0, omega - x,
     ylab = "tpx", xlab = "t")

k <- 0:(omega - x - 1)
dx <- sapply(k, function(t) tpx(t, x, a, omega) * 
               tqx(1, x + t, a, omega)) * radix
plot(x = x + k, y = dx, ylab = "dx", xlab = "x", type = 'o')

k <- 1:(n-1)
ncurtate <- sum(sapply(k, function(t) t * tpx(t, x, a, omega) *
                    tqx(1, x + t, a, omega))) + 
            n * tpx(n, x, a, omega)

ncomplete <- integrate(function(t) tpx(t, x, a, omega), 0, n)$value

medianAge <- head(lifeTable$x[lifeTable$lx < radix/2], 1)
modeAge <- lifeTable$x[which(lifeTable$dx == max(lifeTable$dx))]

sdevTx <- sqrt(2*integrate(function(t) t * tpx(t, x, a, omega), 
                           0, omega - x)$value -
              integrate(function(t) tpx(t, x, a, omega), 
                        0, omega - x)$value^2)

k <- 1:(omega-x-1)
sdevKx <- sqrt(sum(sapply(k, function(t) t^2 * tpx(t, x, a, omega) * 
                            tqx(1, x + t, a, omega))) - 
          sum(sapply(k, function(t) t * tpx(t, x, a, omega) * 
                       tqx(1, x + t, a, omega)))^2)

Rx <- function(x, omega) {
  ifelse(x <= 39, 0.999,
         ifelse(x <= 59, 0.998,
                ifelse(x <= 79, 0.995,
                       ifelse(x <= omega, 0.990, 0))))
}

# s is the number of years to adjust the life table
createAdjustLifeTable <- function(x, radix, a, omega, s) {
  
  k = 0:(omega - x - 1)
  qx <- sapply(k, function(t) tqx(1, x + t, a, omega)) * 
    sapply(x, function(y) Rx(y, omega = omega))^s
  px <- 1 - qx
  
  lx <- radix * c(1, head(cumprod(px), -1))
  dx = c(-diff(lx), tail(lx, 1))
  data.frame(x = x + k, lx = lx, dx = dx, qx = qx, px = px)
}

adjustLifeTable10 <- createAdjustLifeTable(x, radix, a, 
                                           omega, 10)
medianAge10 <- head(adjustLifeTable10$x[adjustLifeTable10$lx < 
                                          radix/2], 1)
modeAge10 <- adjustLifeTable10$x[which(adjustLifeTable10$dx == 
                                         max(adjustLifeTable10$dx))]

adjustLifeTable25 <- createAdjustLifeTable(x, radix, a, 
                                           omega, 25)
medianAge25 <- head(adjustLifeTable25$x[adjustLifeTable25$lx < 
                                          radix/2], 1)
modeAge25 <- adjustLifeTable25$x[which(adjustLifeTable25$dx == 
                                         max(adjustLifeTable25$dx))]