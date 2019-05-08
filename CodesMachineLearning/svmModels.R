library(MASS)
library(e1071)

# Breast cancer data; Diagnosing breast cancer cases as M = malignant or B = benign

data = read.table('breastCancer.txt', sep=',', header=TRUE)
n = dim(data)[1]
p = dim(data)[2]

ind.tr <- sample(n, floor(2*n/3))
ind.te <- setdiff(seq(1, n), ind.tr)

svm.m <- svm(y~., kernel='linear', data=data[ind.tr, ])
pred.class <- predict(svm.m, newdata=data[ind.te, ])
act.class <- data$y[ind.te]
table(act.class, pred.class)
mean(act.class == pred.class)
plot(svm.m, data[ind.te, ], formula= Comp.1~Comp.2)


svm.m <- svm(y~., kernel='polynomial', data=data[ind.tr, ])
pred.class <- predict(svm.m, newdata=data[ind.te, ])
act.class <- data$y[ind.te]
table(act.class, pred.class)
mean(act.class == pred.class)
plot(svm.m, data[ind.te, ], formula= Comp.1~Comp.2)


svm.m <- svm(y~., kernel='radial', data=data[ind.tr, ])
pred.class <- predict(svm.m, newdata=data[ind.te, ])
act.class <- data$y[ind.te]
table(act.class, pred.class)
mean(act.class == pred.class)
plot(svm.m, data[ind.te, ], formula= Comp.1~Comp.2)





# Pima Indian

data = read.table('Pima.txt', header=TRUE)
n = dim(data)[1]
p = dim(data)[2]

ind.tr <- sample(n, floor(2*n/3))
ind.te <- setdiff(seq(1, n), ind.tr)

svm.m <- svm(diabetes ~., kernel='linear', data=data[ind.tr, ])
pred.class <- predict(svm.m, newdata=data[ind.te, ])
act.class <- data$diabetes[ind.te]
table(act.class, pred.class)
mean(act.class == pred.class)

svm.m <- svm(diabetes ~., kernel='radial', data=data[ind.tr, ])
pred.class <- predict(svm.m, newdata=data[ind.te, ])
act.class <- data$diabetes[ind.te]
table(act.class, pred.class)
mean(act.class == pred.class)


# You can also use the ksvm function from kernlab
library(kernlab)
svm.m <- ksvm(diabetes~bmi+glu, type = 'C-svc',  kernel = 'rbfdot', prob.model = TRUE,  data=data)

pdf('PimaSVM.pdf')
plot(svm.m, data=data[, c('bmi', 'glu')])
legend('topleft', legend=c('Yes', 'No'), pch=c(1, 2))
dev.off()




# Letter recognision 

data = read.table('letter.txt', sep=',', header=F)
n = dim(data)[1]
p = dim(data)[2]

ind.tr <- sample(n, floor(2*n/3))
ind.te <- setdiff(seq(1, n), ind.tr)

colnames(data)[1]<-'y'


svm.m <- svm(y~., kernel='linear', data=data[ind.tr, ])
pred.class <- predict(svm.m, newdata=data[ind.te, ])
act.class <- data$y[ind.te]
mean(act.class == pred.class)


svm.m <- svm(y~., kernel='radial', data=data[ind.tr, ])
pred.class <- predict(svm.m, newdata=data[ind.te, ])
act.class <- data$y[ind.te]
mean(act.class == pred.class)
