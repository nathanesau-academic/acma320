# Acma 320 - Spring 2014
# Assignment 7
# Nathan Esau

d <- 2
A <- 0.00022
B <- 2.5e-05
c <- 1.1
omega <- 131
radix <- 1e+05
x <- 40
i <- 0.05
n <- 10
m <- 4

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

v <- function(i, n) {
  (1+i)^(-n)
}

tEx <- function(t, x, s, i) {
  v(i, t) * tpx(t, x, s) 
}

createAnnuityTable <- function(x, radix, A, B, c, omega, d, i, 
                                 moment) {
  
  i <- exp(moment * log(1+i)) - 1
  
  Ax = vector("list", d + 1) # whole life insurance
  ax = vector("list", d + 1) # annuity table
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
  }
  
  ax <- lapply(Ax, function(y) (1 - y)/(i/(1+i)))
  
  # combine lists into data frame (insurance table)
  it = data.frame(x:(omega-d), Ax, ax, x:(omega-d) + d) 
  
  # rename data frame
  if(d>0) { 
    names(it) = c("x", paste0(ifelse(moment == 1, "", moment), "A[x]"), 
                  sapply(1:(d-1), function(d) {
                    paste0(ifelse(moment == 1, "", moment), "A[x]+", d)
                  }), 
                  paste0(ifelse(moment == 1, "", moment), "Ax+", d), 
                  paste0(ifelse(moment == 1, "", moment), "a[x]"), 
                  sapply(1:(d-1), function(d) {
                    paste0(ifelse(moment == 1, "", moment), "a[x]+", d)
                  }), 
                  paste0(ifelse(moment == 1, "", moment), "ax+", d),
                  paste0("x+",d))
  } else {
    names(it) = c("x", "Ax", "ax", "x")
    if(moment > 1) {
      names(it)[2:3] = paste0(moment, ":", names(it)[2:3])
    }
  }
  
  # remove last row
  head(it,-1)
}

annuityTable <- createAnnuityTable(x, radix, A, B, c, 
                                       omega, d, i, 1)

udeferredtqx <- function(u, t, x, s) {
  tpx(u, x, s) - tpx(u + t, x, s)
}


nEx <- function(n, x, s, i, moment) {
  i <- exp(moment * log(1+i)) - 1
  
  tpx(n, x, s) * v(i, n)
}

Ax <- function(x, s, i, moment, n, m, type = 'term') {
  ioriginal <- i
  i <- exp(moment * log(1+i)) - 1
  iupperm <- m*((1+i)^(1/m) - 1) # UDD

  k <- seq(0, pmin(n, omega-x) - 1, 1)
  pk <- udeferredtqx(k, 1, x, s)
  vk <- v(i, k+1)
  EPVterm <- sum(pk * vk) * i/iupperm
  
  if(type == 'endowment') {
    return(EPVterm + nEx(n, x, s, ioriginal, moment))
  } else { # term of or whole life (for whole life use large n)
    return(EPVterm)
  }
}

dupperm <- function(i, m) {
  dannual <- i/(1+i)
  m*(1 - (1-dannual)^(1/m))
}

annx <- function(x, s, i, moment, n, m) {
  (1 - Ax(x, s, i, moment, n, m, type='endowment'))/dupperm(i,m)
}

EPVInsurance1 <- annx(x, 0, i, 1, n, m)
VarInsurance1 <- (Ax(x, 0, i, 2, n, m, type = 'endowment') - 
                  Ax(x, 0, i, 1, n, m, type = 'endowment')^2) * dupperm(i, m)^(-2)

EPVInsurance2 <- nEx(n, x, 0, i, 1) * (annx(x + n, n, i, 1, omega - x - 1, m))
VarInsurance2 <- nEx(n, x, 0, i, 2) * annx(x + n, n, i, 1, omega - x - 1, m)^2 -
                 (nEx(n, x, 0, i, 1) * annx(x + n, n, i, 1, omega - x - 1, m))^2 +
                 nEx(n, x, 0, i, 2) * (Ax(x + n, n, i, 2, omega - x - 1, m) - 
                    Ax(x + n, n, i, 1, omega - x - 1, m)^2) * dupperm(i, m)^(-2)

Percentile <- function(EPV, Var, p, c) {
  qnorm(p) * sqrt(Var) * sqrt(c) + c * EPV
}

p70Insurance1_100 <- Percentile(EPVInsurance1, VarInsurance1, 0.70, 100)
p95Insurance1_100 <- Percentile(EPVInsurance1, VarInsurance1, 0.95, 100)
p70Insurance2_100 <- Percentile(EPVInsurance2, VarInsurance2, 0.70, 100)
p95Insurance2_100 <- Percentile(EPVInsurance2, VarInsurance2, 0.95, 100)

p70Insurance1_10000 <- Percentile(EPVInsurance1, VarInsurance1, 0.70, 10000)
p95Insurance1_10000 <- Percentile(EPVInsurance1, VarInsurance1, 0.95, 10000)
p70Insurance2_10000 <- Percentile(EPVInsurance2, VarInsurance2, 0.70, 10000)
p95Insurance2_10000 <- Percentile(EPVInsurance2, VarInsurance2, 0.95, 10000)