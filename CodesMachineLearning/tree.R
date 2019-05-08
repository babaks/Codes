library(pROC)
library(boot)
library(randomForest)
library(rpart)
library(party)
library(car)

# We use the Mroz dataset for US women's labor-force participation
help(Mroz)
head(Mroz)


###### RPART

# Classification tree to predict the work status
fit.part <- rpart(lfp ~ k5+k618+age+wc+hc+inc, data=Mroz)

plot(fit.part, uniform=TRUE, margin = 0.1)
text(fit.part, use.n=TRUE, cex=.6)

pfit<- prune(fit.part, cp=fit.part$cptable[which.min(fit.part$cptable[,"xerror"]),"CP"])
# plot(pfit, uniform=TRUE, margin = 0.1)
# text(pfit, use.n=TRUE, all=TRUE, cex=.55)

pdf('MrozTree.pdf')
plot(pfit, uniform=TRUE, margin = 0.1)
text(pfit, fancy = TRUE, fheight = 0.3, fwidth = 0.3, use.n=TRUE, all=TRUE, cex=.6)
dev.off()

# Regression tree to predict the income among working women
fit.part <- rpart(lwg ~ k5+k618+age+wc+hc+inc, data=Mroz[Mroz$lfp=='yes', ])

plot(fit.part, uniform=TRUE, margin = 0.1)
text(fit.part, use.n=TRUE, cex=.6)

pfit<- prune(fit.part, cp=fit.part$cptable[which.min(fit.part$cptable[,"xerror"]),"CP"])
plot(pfit, uniform=TRUE, margin = 0.1)
text(pfit, use.n=TRUE, all=TRUE, cex=.55)


###### Random Forest for classification

fit.rf <- randomForest(lfp ~ k5+k618+age+wc+hc+inc, data=Mroz, ntree=1000)

varImpPlot(fit.rf)


# Training and test
n = nrow(Mroz)

ind.tr <- sample(n, floor(2*n/3))
ind.te <- setdiff(seq(1, n), ind.tr)

fit.rf <- randomForest(lfp ~ k5+k618+age+wc+hc+inc, data=Mroz[ind.tr, ], ntree=1000)

y.hat <- predict(fit.rf, newdata=Mroz[ind.te, ], type="response")
tbl<-table(Mroz$lfp[ind.te], y.hat)

acc <- sum(diag(tbl))/sum(tbl)

p <- predict(fit.rf, newdata=Mroz[ind.te, ], type="prob")[, 2]

roc(Mroz$lfp[ind.te], p, plot=TRUE)

auc <- roc(Mroz$lfp[ind.te], p, plot=FALSE)$auc


