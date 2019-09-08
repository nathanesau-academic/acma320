d <- 2
A <- 0.00022
B <- 2.5e-05
c <- 1.1
omega <- 131
radix <- 1e+05
x <- 20
i <- 0.05

tpx <- function(t,x,s) {
  
  f <- function(t,x,s) {
    u = pmax(0,d-s-t)
    j = t <= (d-s)
    
    exp(0.9^(u) * (
      A*(1-j*0.9^(t))*t^(1-j)/(log(0.9)^j - 2*(1 - j)) +
        B*c^(x)*(c^(t) - 0.9^(j*t))/(log(0.9)*j - log(c))
    ))
    
  }
  
  m = pmax(0,d-s)
  f(m,x,s)* f(t-m,x+m,d)
}

createLifeTable <- function(x, radix, A, B, c, omega, d) {
  
  if(d > 0) {
    
    # creates the select life table
    lt = data.frame(
      c(rep(NA, d), x:(omega-d-1)), 
      lapply(0:(d-1), 
            function(y) {
              c(rep(NA, d), 
                tail(tpx(0:(omega-x-1), x, s = d)*radix, -d) / 
                sapply(x:(omega-d-1), 
                       function(x) {
                          tpx(d - y, x + y, y)
                       }
                       ))
            }
      ),
      tpx(0:(omega-x-1), x, s = d)*radix, x:(omega-1)
    )
    
    # renames the select life table
    names(lt) =  c("x", sapply(0:(d-1), function(x) {
                                          paste0("l[x]+",x)
                                        }), 
                    paste0("lx+",d), paste0("x+",d))
  } else { 
    
    # d = 0
    lt = data.frame(x:(omega-d-1), tpx(0:(omega-x-1), x - d, s = d)*radix,
                    x:(omega-1))
    names(lt) = c("x","lx","x")
  }
  
  lt
}

lifeTable <- createLifeTable(x, radix, A, B, c, omega, d)

v <- function(i, n) {
  (1+i)^(-n)
}

tEx <- function(t, x, s, i) {
  v(i, t) * tpx(t, x, s) 
}

createInsuranceTable <- function(x, radix, A, B, c, omega, d, i, n, 
                                 moment) {
  
  i <- exp(moment * log(1+i)) - 1
  
  Ax = vector("list", d + 1) # whole life insurance
  Ex = vector("list", d + 1) # pure endowment insurance
  p = vector("list", d + 1)
  
  # recursive insurance
  recins <- function(p, init = FALSE, Ax = NA) { 
    
    # Calculate A[x] values using recursion
    prev = v(i,1)
    A = prev
    
    # Special case for init
    for(t in (omega-d-1):x) {
      
      prev = ifelse(init, (1-p[t-x+1])*v(i,1) + p[t-x+1]*v(i,1)*prev, 
                    (1-p[t-x+1])*v(i,1) + p[t-x+1]*v(i,1)*Ax[t-x+1])
      A = c(A, prev)
    } 
    
    # A is backwards (since we use a backwards recursion). 
    # Need to reverse A
    rev(A)
  }
  
  # calculate probabilities of survival over one year intervals
  p = lapply(0:d, function(t) {
                    sapply((omega-d-1):x, function(s) {
                                        tpx(1, s + t, t)
                                      }
                    )
                  }
  )
  
  # calculate select insurance formulas recursively
  for(t in d:0) {
    
    if(t == d) {
      Ax[[t+1]] = recins(rev(p[[t+1]]), TRUE) 
    } else {
      Ax[[t+1]] = recins(rev(p[[t+1]]), FALSE, Ax[[t+2]])
    }
    
    Ex[[t+1]] = sapply(x:(omega-d), function(s) { 
      tEx(n,s,t,i)
    }
    )
  }
  
  # combine lists into data frame (insurance table)
  it = data.frame(x:(omega-d), Ax, Ex, x:(omega-d) + d) 
  
  # rename data frame
  if(d>0) { 
    names(it) = c("x", paste0(ifelse(moment == 1, "", moment), "A[x]"), 
      sapply(1:(d-1), function(d) {
        paste0(ifelse(moment == 1, "", moment), "A[x]+", d)
      }), 
      paste0(ifelse(moment == 1, "", moment), "Ax+", d), 
      paste0(ifelse(moment == 1, "", paste0(moment,":")), 
             paste0(n,"E[x]")), 
      sapply(1:(d-1), function(d) {
        paste0(ifelse(moment == 1, "", paste0(moment, ":")), 
               paste0(n,"E[x]+"), d)
      }),
      paste0(ifelse(moment == 1, "", paste0(moment, ":")), 
             paste0(n, "Ex+"), d), paste0("x+",d))
  } else {
    names(it) = c("x", "Ax", paste0(n,"Ex"), "x")
    if(moment > 1) {
      names(it)[2:3] = paste0(moment, ":", names(it)[2:3])
    }
  }
  
  # remove last row
  head(it,-1)
}

insuranceTable <- createInsuranceTable(x, radix, A, B, c, 
                                       omega, d, i, 5, 1)

Ax4 <- createInsuranceTable(x, radix, A, B, c, omega, d, 0.04, 
                            5, 1)$"A[x]"
Ax5 <- createInsuranceTable(x, radix, A, B, c, omega, d, 0.05,
                            5, 1)$"A[x]"
Ax6 <- createInsuranceTable(x, radix, A, B, c, omega, d, 0.06,
                            5, 1)$"A[x]"

plot(x = x:(omega-d-1), Ax4, lty = 1, type = 'l',
    ylab = "A[x]", xlab = "x", ylim = c(min(Ax6), max(Ax4)))
lines(x = x:(omega-d-1), Ax5, lty = 2, type = 'l')
lines(x = x:(omega-d-1), Ax6, lty = 3, type = 'l')

legend('topleft', c("i = 4%", "i = 5%", "i = 6%"),
       lty = c(1, 2, 3))

udeferredtqx <- function(u, t, x, s) {
  tpx(u, x, s) - tpx(u + t, x, s)
}

createProbTable <- function(x, i) {
  
  k <- 0:(omega-x-1)
  pk <- udeferredtqx(k, 1, x, 0)
  Fk <- cumsum(pk)
  vk <- v(i,k)
  pz <- c(0, head(pk, -1))
  Fz <- cumsum(pz)
  
  data.frame(k = k, pk = pk, Fk = Fk,
             vk = vk, pz = pz, Fz = Fz)
}

probTable <- createProbTable(50, i)

Ax <- function(x, s, i, moment) {
  i <- exp(moment * log(1+i)) - 1
  k <- 0:(omega-x-1)
  pk <- udeferredtqx(k, 1, x, s)
  vk <- v(i, k+1)
  sum(pk * vk)
}

ExpectedValueZ <- Ax(50, 0, i, 1)
StandardDeviationZ <- sqrt(Ax(50, 0, i, 2) - Ax(50, 0, i, 1)^2)

Fz <- function(z, x, s, i) {
  k <- floor(uniroot(function(t) v(i,t) - z, 
                     interval = c(0, omega - x - 1))$root)
  1 - sum(udeferredtqx(k:(omega-x-1), 1, x, s))
}

FExpectedValueZ <- Fz(ExpectedValueZ, 50, 0, i)

SimulateK <- function(n, x, omega) {
  k <- 0:(omega-x-1)
  pk <- udeferredtqx(k, 1, x, 0)
  Fk <- cumsum(pk)
  
  u <- runif(n)
  sapply(u, function(y) head(k[Fk > y], 1)) + 1
}

set.seed(2000)

kvalues <- SimulateK(50, 50, omega)
zvalues <- v(i, kvalues)

hist(zvalues, main = "", xlab = "Simulated value of Z")

EstimatedExpectedValueZ <- mean(zvalues)
EstimatedStandardDeviationZ <- sd(zvalues)
EstimatedFExpectedValueZ <- sum(zvalues > ExpectedValueZ) / 50

kvalues <- SimulateK(500, 50, omega)
zvalues <- v(i, kvalues)

hist(zvalues, main = "", xlab = "Simulated value of Z")

EstimatedExpectedValueZ <- mean(zvalues)
EstimatedStandardDeviationZ <- sd(zvalues)
EstimatedFExpectedValueZ <- sum(zvalues > ExpectedValueZ) / 500
