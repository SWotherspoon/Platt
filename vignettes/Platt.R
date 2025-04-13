## -----------------------------------------------------------------------------
library(Platt)
data(PI)
head(PI)

## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(PI,aes(x=I,y=P,group=Depth))+
  geom_point()+
  facet_wrap(~Depth)

## -----------------------------------------------------------------------------
PI.C <- subset(PI,Depth=="C")

## -----------------------------------------------------------------------------
fit <- nls(P~SSPlatt(I,alpha,beta,Pmax,R),data=PI.C)
summary(fit)

## -----------------------------------------------------------------------------
fitted(fit)

## -----------------------------------------------------------------------------
PI.pr <- data.frame(I=seq(0,1700,20))
PI.pr$P <- predict(fit,PI.pr)
head(PI.pr)

## -----------------------------------------------------------------------------
with(as.list(coef(fit)),
  correctedPmax(alpha,beta,Pmax))

## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(PI.C,aes(x=I,y=P))+
  geom_line(data=PI.pr,col="dodgerblue")+
  geom_point()+
  geom_point(aes(x=I,y=fitted(fit)),col="dodgerblue")

## -----------------------------------------------------------------------------
plot(residuals(fit)~fitted(fit))

## -----------------------------------------------------------------------------
by(PI,PI$Depth,function(d) summary(nls(P~SSPlatt(I,alpha,beta,Pmax,R),data=d)))

## -----------------------------------------------------------------------------
library(dplyr)
PI.coef <- PI %>% 
  group_by(Depth) %>% 
  do(as.data.frame(t(coef(nls(P~SSPlatt(I,alpha,beta,Pmax,R),data=.))))) %>%
  as.data.frame

PI.coef

## -----------------------------------------------------------------------------
library(dplyr)
PI.pr <- PI %>% 
  group_by(Depth) %>% 
  do({
    d <- data.frame(I=seq(0,1700,20))
    d$P <- predict(nls(P~SSPlatt(I,alpha,beta,Pmax,R),data=.),d)
    d
  })
ggplot(PI,aes(x=I,y=P,group=Depth))+
  geom_line(data=PI.pr)+
  geom_point()+
  facet_wrap(~Depth)

