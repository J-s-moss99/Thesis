## identifying life-history parameters of the Thornback ray using tagging data

# load packages for plotting and date-time operations
library(lattice)
library(lubridate)
library(broom)
library(tidyverse)

# load data

load("tagging_data.RData")

#################
## 1. Maturity ##
#################

# select the Thornback ray
tagging_data.rjc <- tagging_data[tagging_data$Species=="THR",]
levels(tagging_data.rjc$Maturity)

# keep records that have been assigned a maturity stage
# A: not
# B: developing
# C: mature
# D: reproducing

tagging_data.rjc <- subset(tagging_data.rjc, Maturity %in% c("A","B","C","D"))
tagging_data.rjc <- subset(tagging_data.rjc, Sex %in% c("M","F"))

# define binary variable: mature = 1; immature = 0
tagging_data.rjc$mature <- 1
tagging_data.rjc$mature[tagging_data.rjc$Maturity %in% c("A","B")] <- 0 

# some data exploration
xyplot(mature~Length_cm|Sex,data=tagging_data.rjc)


# fit models to data
# binomial data

mod.maturity.L   <- glm(mature ~ Length_cm, family=binomial(link="logit"),data=tagging_data.rjc)

# quick diagnostics
# check influential values
plot(mod.maturity.L, which = 4, id.n = 3)
# Extract model results
model.data <- augment(mod.maturity.L) %>% mutate(index = 1:n()) 
# display the data (top 3) with largest Cooks distance
model.data %>% top_n(3, .cooksd)

ggplot(model.data, aes(index, .std.resid)) + 
  geom_point(aes(color = mature), alpha = .5) +
  theme_bw()


mod.maturity.S   <- glm(mature ~ Length_cm + as.character(Sex), family=binomial(link="logit"),data=tagging_data.rjc)
mod.maturity.SL  <- glm(mature ~ Length_cm * as.character(Sex), family=binomial(link="logit"),data=tagging_data.rjc)

### Model comparison

# compare models based on AIC
AIC(mod.maturity.L);AIC(mod.maturity.S);AIC(mod.maturity.SL) # second model (with length and sex has lowest AIC)

# other comparison could be to use a subset of the data for training, and a subset for validation
#----#

# compare predictions with observed values, in table

tagging_data.rjc$pred.L <- as.numeric(predict(mod.maturity.L,newdata =tagging_data.rjc,type="response",se.fit = T)$fit >0.5)
sum(table(tagging_data.rjc$pred.L,tagging_data.rjc$mature,tagging_data.rjc$Sex))

tagging_data.rjc$pred.S <- as.numeric(predict(mod.maturity.S,newdata = tagging_data.rjc,type="response",se.fit = T)$fit >0.5)
sum(table(tagging_data.rjc$pred.S,tagging_data.rjc$mature,tagging_data.rjc$Sex))

tagging_data.rjc$pred.SL <- as.numeric(predict(mod.maturity.SL,newdata = tagging_data.rjc,type="response",se.fit = T)$fit >0.5)
sum(table(tagging_data.rjc$pred.SL,tagging_data.rjc$mature,tagging_data.rjc$Sex))

# or create a ROC curve to compare the area under the curve
library(pROC)

roccurve.L  <- roc(tagging_data.rjc$mature ~ tagging_data.rjc$pred.L)
roccurve.S  <- roc(tagging_data.rjc$mature ~ tagging_data.rjc$pred.S)
roccurve.SL <- roc(tagging_data.rjc$mature ~ tagging_data.rjc$pred.SL)

auc(roccurve.L)
auc(roccurve.S)    # largest area with fewest parameters
auc(roccurve.SL)

plot(roccurve.L)
plot(roccurve.S,add=T,col="red")
plot(roccurve.SL,add=T,col="blue")

###### plot Pr(maturity) over length

res.pred <- predict(mod.maturity.S,
                    newdata = data.frame(Length_cm=rep(seq(0,120,0.1),2),
                                         Sex=c(rep("M",length(seq(0,120,0.1))),rep("F",length(seq(0,120,0.1))))),
                    type="response",se.fit = T)

# plot and add lines
par(mfrow=c(2,1))
# males
plot(seq(0,120,0.1),res.pred$fit[1:length(seq(0,120,0.1))],type="l",xlab="length (cm)",ylab="Pr(Mature)",main="males")
polygon(c(seq(0,120,0.1),rev(seq(0,120,0.1))),
        c(res.pred$fit[1:length(seq(0,120,0.1))] + 1.96*res.pred$se.fit[1:length(seq(0,120,0.1))] , 
          rev(res.pred$fit[1:length(seq(0,120,0.1))] - 1.96 *res.pred$se.fit[1:length(seq(0,120,0.1))])),
        col="grey60")
lines(seq(0,120,0.1),res.pred$fit[1:length(seq(0,120,0.1))],type="l")
# length at 50% maturity
M50 <- which(abs(res.pred$fit[1:length(seq(0,120,0.1))]-0.5) == min(abs(res.pred$fit[1:length(seq(0,120,0.1))]-0.5)))/10
abline(v=M50,col="red",lty="dotted")
text(M50-20,0.9,
     paste0("M 50% = ",M50),col="red",cex=0.8)

# females
plot(seq(0,120,0.1),res.pred$fit[(length(seq(0,120,0.1))+1):length(res.pred$fit)],type="l",xlab="length (cm)",ylab="Pr(Mature)",main="females")
polygon(c(seq(0,120,0.1),rev(seq(0,120,0.1))),
        c(res.pred$fit[(length(seq(0,120,0.1))+1):length(res.pred$fit)] + 1.96*res.pred$se.fit[(length(seq(0,120,0.1))+1):length(res.pred$fit)] , 
          rev(res.pred$fit[(length(seq(0,120,0.1))+1):length(res.pred$fit)] - 1.96 *res.pred$se.fit[(length(seq(0,120,0.1))+1):length(res.pred$fit)])),
        col="grey60")
lines(seq(0,120,0.1),res.pred$fit[(length(seq(0,120,0.1))+1):length(res.pred$fit)],type="l")
M50 <- which(abs(res.pred$fit[(length(seq(0,120,0.1))+1):length(res.pred$fit)]-0.5) == min(abs(res.pred$fit[(length(seq(0,120,0.1))+1):length(res.pred$fit)]-0.5)))/10
abline(v=M50,col="red",lty="dotted")
text(M50-20,0.9,
     paste0("M 50% = ",M50),col="red",cex=0.8)

