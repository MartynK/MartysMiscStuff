#######
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

#######
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

#######
########## eloszlás

set.seed(1)
ex.hist2 <- ex_dist1(50000)

par(mfcol=c(2,1))
hist(ex.hist2,xlab="",main="~Sűrűség 50k mérés alapján",
     cex.axis=1,cex.lab=1,freq=F,breaks=160)

plot(ecdf(ex.hist2),xlim=c(100,150),
     xlab="",main="(empírikus) eloszlásfüggvény")

par(mfcol=c(1,1))

#######
####### binomial

ex.binom <- rep(0,6)
for (i in 1:6) {
  ex.binom[i] <- dbinom(i-1,5,0.5) 
}

barplot(ex.binom,names.arg=seq(0,5),xlab="Fejek száma",ylab="Valószínűség",
        main="Fejek számának valószínűsége 5 pénzfeldobás után")

#######
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

#######
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

#######
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

#######
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

#######
#### Powsampsize2

curve(dnorm(x, mean=0, sd=3),from=-20,to=30)
curve(dnorm(x, mean=6, sd=3),add=T)

#####
# "szakaszos" ábra két mintára + egy minta jellemzésére

set.seed(123456)
no_points <- 25
a <- sort(rnorm(no_points,-2))
b <- sort(rnorm(no_points,2))


## define plot data
xlim <- c(-5,5);
ylim <- c(0,1.3);
py <- rep(0,length(a));
lx.buf <- 5;
lx <- seq(xlim[1]+lx.buf,xlim[2]-lx.buf,len=length(py));
ly <- 20;

## create basic plot outline
par(mfrow=c(4,1))

par(xaxs='i',yaxs='i',mar=c(5,1,1,1));
plot(NA,xlim=xlim,ylim=ylim,axes=F,ann=F);
axis(1);

## plot elements
points(a,py+0.1,pch=16,xpd=NA,col='blue')
points(b,py+0.2,pch=16,xpd=NA,col='red')

#means
points(mean(a),py[1]+0.05,pch=15,xpd=NA,cex=1.5)
points(mean(b),py[1]+0.25,pch=15,xpd=NA,cex=1.5)

# Single sd()
segments(mean(a),py[1]+0.05,mean(a)-sd(a))
segments(mean(a),py[1]+0.05,mean(a)+sd(a))
segments(mean(b),py[1]+0.25,mean(b)-sd(b))
segments(mean(b),py[1]+0.25,mean(b)+sd(b))

#95% CI
segments(mean(a),py[1]+0.05,mean(a)-(sd(a)*1.96),lty=2)
segments(mean(a),py[1]+0.05,mean(a)+(sd(a)*1.96),lty=2)
segments(mean(b),py[1]+0.25,mean(b)-(sd(a)*1.96),lty=2)
segments(mean(b),py[1]+0.25,mean(b)+(sd(a)*1.96),lty=2)

#Boxplots
boxplot(a,horizontal=T,add=T,border = T,at=0.5)
boxplot(b,horizontal=T,add=T,border = T,at=0.5)
# 
# points(a,py+0.1,pch=16,xpd=NA,col='blue')


sum_a <- summary(a)
(sum_a)

#par(mfrow=c(2,1))
hist(a,xlim=c(-4,1),freq=F)
lines(density(a))
plot(ecdf(a),xlim=c(-4,1))

qqnorm(a)
qqline(a)

par(mfrow=c(1,1))


#####
# Példa átlag szórására

par(mfcol=c(1,1))
set.seed(123456)
no_points <- 25
xlim <- c(-5,3);
ylim <- c(0,1.3);


par(xaxs='i',yaxs='i',mar=c(5,1,1,1));
plot(NA,xlim=xlim,ylim=ylim,axes=F,ann=F);
axis(1);

no_reps <- 1000

means <- rep(NA,no_reps)
for (i in 1:no_reps) {
  
  py <- rep(0+i*(1/no_reps),length(a));
  
  a <- sort(rnorm(no_points,-2))
  
  ## plot elements
  points(a,py+0.1,pch=16,xpd=NA,col='blue')
  
  #means
  points(mean(a),py[1]+0.1,pch=15,xpd=NA,cex=1.5)
  
  means[i] <- mean(a)  
}

hist(means,col="gray")

par(mfcol=c(1,1))

##### 
# Példa "előtte/utána" összehasonlításnál a véletlen szerepére

xlim <- c(-10,20);
ylim <- c(0,1.3);
py <- rep(0,length(a));
lx.buf <- 5;
lx <- seq(xlim[1]+lx.buf,xlim[2]-lx.buf,len=length(py));
ly <- 20;

## create basic plot outline

par(xaxs='i',yaxs='i',mar=c(12,1,1,1));
plot(NA,xlim=xlim,ylim=ylim,axes=F,ann=F);
axis(1);

arrows(x0=0,y0=.5,x1=0,y1=0,lwd=2)

arrows(x0=5,y0=1.5,x1=5,y1=0,lwd=2,col="red")
sd_s  <- 7
n_exa <- 34
a     <- rnorm(500000,0,sd_s)
boxplot(a,horizontal=T,add=T,border = T,at=1,outline=F,range=0.5)

segments(x0=(0-sd_s/sqrt(n_exa))*2,x1=(0+sd_s/sqrt(n_exa))*2,
         y0=0.6,y1=0.6,
         col="darkblue",lwd=4)

segments(x0=5-(sd_s/sqrt(n_exa))*2,x1=5+(sd_s/sqrt(n_exa))*2,
         y0=0.65,y1=0.65,
         col="red",lwd=4)


set.seed(12345)
n_rep   <- 5000
n_exa   <- 43
outcome <- rep(F,n_rep)
for (i in 1:n_rep) {
  samp_a     <- rnorm(n_exa,0,sd_s)
  samp_b     <- rnorm(n_exa,5,sd_s)  
  outcome[i] <- ifelse(t.test(samp_a,samp_b)$p.value<0.05,T,F)
}
sum(outcome)/n_rep

##### 
# példa lm-re, baseline korrekcióra

set.seed(12345)
gen_data <- function(n) {
  x1 <- rnorm(n,100,12)
  x2 <- runif(n,160,190)
  x3 <- rbinom(n,1,0.5)
  y  <- rep(NA,n)
  for (i in 1:n) {
    y_a  <- x1[i]+2*(x2[i]-180)-10*x3[i]
    y[i] <- rnorm(1,y_a,5)
  }
  out <- data.frame(data=y)
  colnames(out)     <- "y"
  out$bp1           <- x1
  out$height        <- x2
  out$treat         <- as.factor(x3)
  return(out)
}

data <- gen_data(300)
pairs(data[,1:3])
plot(data$y,data$bp1)
mod1 <- lm(y~treat,data=data)
summary(mod1)
plot(effects::predictorEffect("treat",mod1))
#plot(mod1)

mod2 <- lm(y~treat+bp1+height,data=data)
summary(mod2)
anova(mod2)
plot(effects::predictorEffect("treat",mod=mod2))
plot(effects::predictorEffect("height",mod=mod2))
plot(effects::predictorEffect("bp1",mod=mod2))
#plot(mod2)

##### 
# Példa torzítatlan szórás hatására
var_samp <- function(input) {
  out <- 0
  ave <- mean(input)
  for (i in 1:length(input)) {
    out <- out +  (input[i]-ave)^2
    
  }
  return(out/length(input))
}

set.seed(12345)
a <- rnorm(8,0,5)
data <- a
sd(data)
sqrt(var_samp(data))

par(xaxs='i',yaxs='i',mar=c(5,1,1,1));
plot(NA,xlim=xlim,ylim=ylim,axes=F,ann=F);
axis(1);

## plot elements
points(a,y=rep(0.1,8),pch=16,xpd=NA,col='blue')

#mean
points(mean(a),py[1]+0.2,pch=15,xpd=NA,cex=1.5)

#95% CI
segments(mean(a),py[1]+0.2,mean(a)-(sd(a)*1.96),lty=2)
segments(mean(a),py[1]+0.2,mean(a)+(sd(a)*1.96),lty=2)

#Biased 95%CI
points(mean(a),py[1]+0.3,pch=15,xpd=NA,cex=1.5,col="red")
segments(mean(a),py[1]+0.3,mean(a)-(sqrt(var_samp(data))*1.96),lty=2,col="red")
segments(mean(a),py[1]+0.3,mean(a)+(sqrt(var_samp(data))*1.96),lty=2,col="red")

par(mfcol=c(1,1))


#####
# Példa normális adatok hisztogramjainak "recegésére"
a <- rnorm(10000,0,1)
#hist(a,probability = T,type="n")
plot(1,1,type="n",xlim=c(-4,4),ylim=c(0,0.5))
lines(density(a),lwd=4)   


##### 
# Becsapós példa: melyik minta különbözik a többitől?
# Válasz: egyik sem

par(mfcol=c(1,1))
set.seed(123456)
no_points <- 8
xlim <- c(80,120);
ylim <- c(0,1.3);


par(xaxs='i',yaxs='i',mar=c(5,1,1,1));
plot(NA,xlim=xlim,ylim=ylim,axes=F,ann=F);
axis(1);

no_reps <- 4

means <- rep(NA,no_reps)
for (i in 1:no_reps) {
  
  py <- rep(0+i*(1/no_reps),length(a));
  
  a <- sort(rnorm(no_points,100,5))
  
  ## plot elements
  points(a,py+0.1,pch=16,xpd=NA,col='blue')
  
}

#####
# "Kontroll limitek" - ismert szórás esetén +-3 szigmán kívül 
# (kb.) 1/1M az esélye a véletlen által okozott előfordulásnak

set.seed(123456)
a <- rnorm(1000,100,5)
plot(a,ylim=c(50,150),xlab="")
abline(h=90,col="blue")
abline(h=110,col="blue")

abline(h=70,col="red")
abline(h=130,col="red") # 6 sigma, 3.4 DPMO

points(420,140, col="red",pch=19)




