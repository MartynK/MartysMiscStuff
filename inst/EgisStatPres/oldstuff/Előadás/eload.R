jegyek <- c(
2,
5,
5,
3,
5,
5,
5,
2,
3,
4,
3,
4,
1,
3,
2,
4)

hist(jegyek)


##### datatable példák
set.seed(1)
ex.wide              <- data.frame(data=matrix(data=NA,nrow=10,ncol=10))
ex.wide[,1]          <- seq(1,10)
ex.long              <- data.frame(data=matrix(data=NA,nrow=90,ncol=4))
ex.long[,1]          <- seq(1,90)
cnam                 <- rep("",10)

for(i in 1:9) {
 cnam[i+1]                    <- paste(i,".obs",sep="")
 dat.act                      <- round(rnorm(10,100-i*5,5+i),0)
 ex.wide[,i+1]                <- dat.act
 ex.long[((i-1)*10+1):(i*10),2] <- i
 ex.long[((i-1)*10+1):(i*10),3] <- dat.act
 ex.long[((i-1)*10+1):(i*10),4] <- seq(1,10)
}
cnam[1]           <- "No.rat"
colnames(ex.wide) <- cnam
colnames(ex.long) <- c("No.Obs","Month","Obs","Rat.no")

###### hist.

ex_dist1 <- function(no) {

 frac1n <- floor(no*0.5)
 frac2n <- floor(no*0.25)
 frac3n <- no-(frac1n+frac2n)
 
 frac1 <- runif(frac1n,100,150)
 frac2 <- exp(rnorm(frac2n,0.1,0.2))*5 + 100
 frac3 <- abs(rnorm(frac3n,130,5))
 
 return(c(frac1,frac2,frac3))
}

set.seed(1)
ex.hist <- ex_dist1(200)

par(mfcol=c(1,1))
hist(ex.hist,xlab="",ylab="Gyakoriság",
     main="hisztogram 200 mérés alapján",cex.axis=1.5,cex.lab=1.5)

arrows(x0=mean(ex.hist),  #Átlag
       y0=30,
       y1=0,
       col="red",
       lwd=7
       )
       
arrows(x0=median(ex.hist), #Medián
       y0=50,
       y1=0,
       length=0.1,
       col="cyan",
       lwd=4
)

arrows(x0=mean(ex.hist)-sd(ex.hist), #SD
       y0=50,
       y1=0,
       length=0.1,
       col="gray",
       lwd=4
)


arrows(x0=quantile(ex.hist,probs=.9), #90th percentile
       y0=40,
       y1=0,
       length=0.2,
       col="blue",
       lwd=5
)

#boxplot(ex.hist,horizontal=T,ylim=c(100,150))

summary(ex.hist)

########## eloszlás

set.seed(1)
ex.hist2 <- ex_dist1(50000)

par(mfcol=c(2,1))
hist(ex.hist2,xlab="",main="~Sűrűség 50k mérés alapján",
     cex.axis=1,cex.lab=1,freq=F,breaks=160)

plot(ecdf(ex.hist2),xlim=c(100,150),
     xlab="",main="(empírikus) eloszlásfüggvény")

par(mfcol=c(1,1))

####### binomial

ex.binom <- rep(0,6)
for (i in 1:6) {
 ex.binom[i] <- dbinom(i-1,5,0.5) 
}

barplot(ex.binom,names.arg=seq(0,5),xlab="Fejek száma",ylab="Valószínűség",
        main="Fejek számának valószínűsége 5 pénzfeldobás után")

###### Stat.exa
set.seed(1)
ex.stat <- rbinom(n=1,size=30,prob=0.55) #18
ll      <- matrix(data=NA,nrow=100,ncol=4)
ll[,1]  <- seq(0,1,length.out=100)
for (i in 1:1000) {
 ll[i,2] <- dbinom(x=18,size=30,prob=ll[i,1])
 ll[i,3] <- log(ll[i,2])
 if (i>1) {ll[i,4] <- ll[i,3]-ll[i-1,3]}
}
colnames(ll) <- c("p","Likelihood","log-likelihood")
par(mfcol=c(3,1))
plot(ll[,1],ll[,2],xlab="feltételezett 'p'",main="Binom modellek likelihoodja a kísérlet függvényében")
plot(ll[,1],ll[,3],xlab="feltételezett 'p'",main="log-likelihood függvény")
plot(ll[,1],ll[,4],xlab="feltételezett 'p'",main="score függvény")
abline(h=0,col="red")

par(mfcol=c(1,1))

cis <- binom::binom.confint(18,30)
plot(x=c(0.4,0.8),y=c(-50,70),type="n",yaxt="n",ylab="",xlab="p",
     main="95% CI-k 11 különböző módszerrel (18 fej 30-ból)")
for (i in 1:11) {
 segments(cis[i,5],0-50+10*i,cis[i,6],0-50+10*i,col="blue",lwd=3) 
}
abline(v=.55,col="red",lwd=2)
abline(v=18/30,col="gray",lty=2,lwd=3)

par(mfcol=c(1,1))
n        <- 1000
ex.stat2 <- matrix(data=c(rbinom(n,30,0.55),rep(NA,2*n)),ncol=3,byrow=F)
ex.stat3 <- rep(0,10000)
plot(x=c(0,1),y=c(-100,100),type="n",yaxt="n",
     main="Ismételt kísérletek, konfidencia intervallumokkal és >95% lefedettségi szint ")
for (i in 1:n) {
 act             <- binom.test(ex.stat2[i,1],30)$conf.int
 ex.stat2[i,2:3] <- act
 segments(act[1],0-100+i*2,act[2],0-100+i*2)
 act             <- round(act*10000)
 ex.stat3[act[1]:act[2]] <- ex.stat3[act[1]:act[2]] + 1
}

abline(v=.5060,col="red")
abline(v=.5940,col="red")

##### teszt stat

out <- c()
for (i in 0:30) {
 outc <- c(rep(1,i),rep(0,30-i)) # Unsorted!
 stat <- t.test(outc,mu=0.5)
 out  <- c(out,stat$statistic) 
  
}
out <- as.data.frame(out)
colnames(out) <- "t_stat"
out$index <- seq(0,1,by=1/29)
plot(x=out$index,y=out$t_stat,xlab="Megfigyelt fejek aránya",ylab="t statisztika",pch=19,cex=1.5)

##### P values

binom.test(18,30,p=0.5)

out <- c()
for (i in 0:30) {
  stat <- binom.test(i,30,p=0.5)
  out  <- c(out,stat$p.value) 
  
}
out <- as.data.frame(out)
colnames(out) <- "p_value"
out$index <- seq(0,30)
plot(x=out$index,y=out$p_value,xlab="Megfigyelt fejek száma 30 dobásból",ylab="p-érték",pch=19,cex=1.5)

### Pow&Sample size

par(mfrow=c(2,1))

out <- c()
for (i in 0:30) {
  stat <- binom.test(i,30,p=0.5)
  out  <- c(out,stat$p.value) 
  
}
out <- as.data.frame(out)
colnames(out) <- "p_value"
out$index <- seq(0,30)
plot(x=out$index/30,y=out$p_value,xlab="Megfigyelt fejek száma 30 dobásból",ylab="p-érték",pch=1,cex=1.5)
segments(x0=out$index[which(out$p_value>0.05)[1]]/30,y0=0.2,
         x1=1-out$index[which(out$p_value>0.05)[1]]/30,
         col="red",lwd=3)

out <- c()
for (i in 0:90) {
  stat <- binom.test(i,90,p=0.5)
  out  <- c(out,stat$p.value) 
  
}
out <- as.data.frame(out)
colnames(out) <- "p_value"
out$index <- seq(0,90)
plot(x=out$index/90,y=out$p_value,xlab="Megfigyelt fejek száma 90 dobásból",ylab="p-érték",pch=1,cex=1.5)
segments(x0=out$index[which(out$p_value>0.05)[1]]/90,y0=0.2,
         x1=1-out$index[which(out$p_value>0.05)[1]]/90,
         col="red",lwd=3)


par(mfrow=c(1,1))

#### Powsampsize2

curve(dnorm(x, mean=0, sd=3),from=-20,to=30)
curve(dnorm(x, mean=6, sd=3),add=T)






