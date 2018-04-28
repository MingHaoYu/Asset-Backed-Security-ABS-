#title : 405-Computational Methods HW9
#author: Ming-Hao Yu
#date  : 2018-03-12
WAC <- 0.08
T <- 30
Loan <- 100000
r0 <- 0.078
kappa <- 0.6
rbar <- 0.08
sigma <- 0.12

Problem1 <- function(Loan, r0, T, kappa, rbar, sigma, WAC, seed=10) {
    simulation <- 1000
    dt <- 1/1200
    steps <- T/dt
    r <- WAC/12
    R <- r*12
    N <- 360
    set.seed(seed)
    
    dw <- sqrt(dt)*rnorm(simulation*steps)
    dwTable <- matrix(nrow=simulation, ncol=steps, dw, byrow=T)
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    PVTable <- matrix(nrow=simulation, ncol=N+1, Loan)
    cashFlow <- matrix(nrow=simulation, ncol=N)
    dtTable <- matrix(nrow=simulation, ncol=N)
    h1 = sqrt(kappa^2+2*sigma^2)
    h2 = (kappa+h1)/2
    h3 = (2*kappa*rbar)/sigma^2
    B = (exp(h1*10)-1)/(h2*(exp(h1*10)-1)+h1)
    A = (h1*exp(h2*10)/(h2*(exp(h1*10)-1)+h1))^h3
    
    for(i in 1:steps) {
        #Use Full Truncation to deal with negative r
        rTable[, i+1] <- rTable[, i] + kappa*(rbar-pmax(rTable[, i], 0))*dt + sigma*sqrt(pmax(rTable[, i], 0))*dwTable[, i]
    }
    rCumTable <- matrix(nrow=simulation, ncol=steps)
    for(i in 1:simulation) {
        rCumTable[i, ] <- exp(-cumsum(rTable[i, 1:steps])*dt)
    }
    
    r10 <- RI <- BU <- IP <- CPR <- TPP <- matrix(nrow=simulation, ncol=360)
    SG <- pmin(1, seq(1,360)/30)
    SY <- rep(c(0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.1, 1.18, 1.22, 1.23, 0.98), T) #monthly calculated
    monthly <- seq(T/dt/360, steps, by=T/dt/360)
    dtTable <- rCumTable[, monthly]
    for(i in 1:360) {
        #r10[, i] <- -rowSums(rTable[, (i):(i+120)])*dt/10
        IP[,i] = r*PVTable[,i]
        P = A*exp(-B*rTable[, monthly[i]])
        r10[, i] <- log(1/P)/10
        RI[, i] <- 0.28+0.14*atan(-8.57+430*(R-r10[, i]))
        BU[, i] <- 0.3+0.7*PVTable[, i]/Loan
        CPR[, i] <- RI[, i]*BU[, i]*SG[i]*SY[i]
        TPP[, i] <- PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1)+(PVTable[, i]-PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1))*(1-(1-CPR[, i])^(1/12))
        PVTable[, i+1] <- PVTable[, i] - TPP[, i]
        cashFlow[, i] <- PVTable[, i]*r/(1-(1+r)^(-N+i-1))+(PVTable[, i]-PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1))*(1-(1-CPR[, i])^(1/12))
    }
    price <- mean(rowSums(dtTable*cashFlow))
    IO <- mean(rowSums(dtTable*IP))
    PO <- mean(rowSums(dtTable*TPP))
    return(list(price=price, IO=IO, PO=PO))
}

Problem2 <- function(Loan, r0, T, kappa, rbar, sigma, WAC, seed=10) {
    simulation <- 1000
    dt <- 1/1200
    steps <- T/dt
    r <- WAC/12
    R <- r*12
    N <- 360
    set.seed(seed)
    
    dw <- sqrt(dt)*rnorm(simulation*steps)
    dwTable <- matrix(nrow=simulation, ncol=steps, dw, byrow=T)
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    PVTable <- matrix(nrow=simulation, ncol=N+1, Loan)
    cashFlow <- matrix(nrow=simulation, ncol=N)
    dtTable <- matrix(nrow=simulation, ncol=N)
    rCumTable <- matrix(nrow=simulation, ncol=steps)
    
    for(i in 1:steps) {
        #Use Full Truncation to deal with negative r
        rTable[, i+1] <- rTable[, i] + kappa*(rbar-pmax(rTable[, i], 0))*dt + sigma*sqrt(pmax(rTable[, i], 0))*dwTable[, i]
    }
    for(i in 1:simulation) {
        rCumTable[i, ] <- exp(-cumsum(rTable[i, 1:steps])*dt)
    }
    monthly <- seq(T/dt/360, steps, by=T/dt/360)
    dtTable <- rCumTable[, monthly]
    
    CPR <- pmin(0.2*seq(1,N), 6)/100
    TPP <- matrix(nrow=simulation, ncol=N)
    for(i in 1:N) {
        TPP[, i] <- PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1)+(PVTable[, i]-PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1))*(1-(1-CPR[i])^(1/12))
        PVTable[, i+1] <- PVTable[, i] - TPP[, i]
        cashFlow[, i] <- PVTable[, i]*r/(1-(1+r)^(-N+i-1))+(PVTable[, i]-PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1))*(1-(1-CPR[i])^(1/12))
    }
    price <- mean(rowSums(dtTable*cashFlow))
    return(price)
}

Problem3 <- function(Loan, r0, OAS, T, kappa, rbar, sigma, WAC, seed=10) {
    simulation <- 1000
    dt <- 1/1200
    steps <- T/dt
    r <- WAC/12
    R <- r*12
    N <- 360
    set.seed(seed)
    
    dw <- sqrt(dt)*rnorm(simulation*steps)
    dwTable <- matrix(nrow=simulation, ncol=steps, dw, byrow=T)
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    PVTable <- matrix(nrow=simulation, ncol=N+1, Loan)
    cashFlow <- matrix(nrow=simulation, ncol=N)
    dtTable <- matrix(nrow=simulation, ncol=N)
    h1 = sqrt(kappa^2+2*sigma^2)
    h2 = (kappa+h1)/2
    h3 = (2*kappa*rbar)/sigma^2
    B = (exp(h1*10)-1)/(h2*(exp(h1*10)-1)+h1)
    A = (h1*exp(h2*10)/(h2*(exp(h1*10)-1)+h1))^h3
    
    for(i in 1:steps) {
        #Use Full Truncation to deal with negative r
        rTable[, i+1] <- rTable[, i] + kappa*(rbar-pmax(rTable[, i], 0))*dt + sigma*sqrt(pmax(rTable[, i], 0))*dwTable[, i]
    }
    #rTable <- rTable + OAS
    
    rCumTable <- matrix(nrow=simulation, ncol=steps)
    for(i in 1:simulation) {
        rCumTable[i, ] <- exp(-cumsum(rTable[i, 1:steps]+OAS)*dt)
    }
    
    r10 <- RI <- BU <- CPR <- TPP <- matrix(nrow=simulation, ncol=360)
    SG <- pmin(1, seq(1,360)/30)
    SY <- rep(c(0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.1, 1.18, 1.22, 1.23, 0.98), T) #monthly calculated
    monthly <- seq(T/dt/360, steps, by=T/dt/360)
    dtTable <- rCumTable[, monthly]
    for(i in 1:360) {
        #r10[, i] <- -rowSums(rTable[, (i):(i+120)])*dt/10
        P = A*exp(-B*rTable[, monthly[i]])
        r10[, i] <- log(1/P)/10
        RI[, i] <- 0.28+0.14*atan(-8.57+430*(R-r10[, i]))
        BU[, i] <- 0.3+0.7*PVTable[, i]/Loan
        CPR[, i] <- RI[, i]*BU[, i]*SG[i]*SY[i]
        TPP[, i] <- PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1)+(PVTable[, i]-PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1))*(1-(1-CPR[, i])^(1/12))
        PVTable[, i+1] <- PVTable[, i] - TPP[, i]
        cashFlow[, i] <- PVTable[, i]*r/(1-(1+r)^(-N+i-1))+(PVTable[, i]-PVTable[, i]*r*(1/(1-(1+r)^(-N+i-1))-1))*(1-(1-CPR[, i])^(1/12))
    }
    price <- mean(rowSums(dtTable*cashFlow))
    return(price)
}
#Problem1 (a)
pricea <- Problem1(Loan, r0, T, kappa, rbar, sigma, WAC, 10)$price

#Problem1 (b)
kappaSeq <- seq(0.3, 0.9, 0.1)
priceb <- 0
for(i in 1:length(kappaSeq)) {
    priceb[i] <- Problem1(Loan, r0, T, kappaSeq[i], rbar, sigma, WAC)$price
}
plot(y=priceb, x=kappaSeq, ylab="Price", xlab="kappa", main="Problem1 (b)", type="l")
names(priceb) <- kappaSeq

#Problem1 (c)
rbarSeq <- seq(0.03, 0.09, 0.01)
pricec <- 0
for(i in 1:length(rbarSeq)) {
    pricec[i] <- Problem1(Loan, r0, T, kappa, rbarSeq[i], sigma, WAC)$price
}
plot(y=pricec, x=rbarSeq, ylab="Price", xlab="rbar", main="Problem1 (c)", type="l")
names(pricec) <- rbarSeq

#Problem2 (a)
price2a <- Problem2(Loan, r0, T, kappa, rbar, sigma, WAC)

#Problem2 (b)
kappaSeq <- seq(0.3, 0.9, 0.1)
price2b <- mapply(FUN=Problem2, 
                  rep(Loan, length(kappaSeq)), 
                  rep(r0, length(kappaSeq)),
                  rep(T, length(kappaSeq)),
                  kappaSeq,
                  rep(rbar, length(kappaSeq)),
                  rep(sigma, length(kappaSeq)),
                  rep(WAC, length(kappaSeq)))
plot(y=price2b, x=kappaSeq, ylab="Price", xlab="kappa", main="Problem2 (b)", type="l")
names(price2b) <- kappaSeq

#Problem3
count <- 0
increment <- -0.001
OAS <- -0.012
while(count < 4) {
    for(i in 1:200) {
        p <- Problem3(Loan, r0, OAS, T, kappa, rbar, sigma, WAC)
        if(p>110000) {
            OAS <- OAS - increment*0.9
            count <- count + 1
            increment <- increment*0.1
            break
        }
        print(c(OAS,p, count))
        OAS <- OAS+increment
    }
}
OAS

#Problem4
y <- 0.0005
P0 <- Problem3(Loan, r0, OAS, T, kappa, rbar, sigma, WAC)
P1 <- Problem3(Loan, r0, OAS-y, T, kappa, rbar, sigma, WAC)
P2 <- Problem3(Loan, r0, OAS+y, T, kappa, rbar, sigma, WAC)

OASDuration <- (P1-P2)/(2*y*P0)
OASConvexity <- (P2-2*P0+P1)/(2*y^2*P0)

#Problem5
ans5 <- matrix(ncol=length(rbarSeq), nrow=2)
for(i in 1:length(rbarSeq)) {
    a <- Problem1(Loan, r0, T, kappa, rbarSeq[i], sigma, WAC)
    ans5[1, i] <- a$IO
    ans5[2, i] <- a$PO
}
colnames(ans5) <- rbarSeq
rownames(ans5) <- c("IO", "PO")
ans5
ylim <- c(min(ans5), max(ans5))
plot(y=ans5[1, ], x=rbarSeq, main="Problem 5", xlab="rbar", ylab="Value", type="l", col="blue", ylim=ylim)
points(y=ans5[2, ], x=rbarSeq, type="l", col="red")
legend("topright", c("IO", "PO"), col=c("blue", "red"), lwd=2.5, cex=0.6)

