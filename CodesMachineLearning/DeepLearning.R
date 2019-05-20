require(pROC)
require(h2o)
require(data.table)
h2o.init(nthreads = -1)

# Reading in the breast cancer data. The outcome is a binary variable: M = malignant, B = benign
data = read.table('breastCancer.txt', sep=',', header=TRUE)
n = dim(data)[1]
p = dim(data)[2]

# Divide the data into training and test
ind.tr <- sample(n, floor(2*n/3))
ind.te <- setdiff(seq(1, n), ind.tr)

TR <- data[ind.tr, ]
TE <- data[ind.te, ]

# Defining the datasets at h2o objects
train <- as.h2o(TR)
test <- as.h2o(TE)

x.indep <- 1:10 # Indeces of the predictors
y.dep <- 11 # Index of the outcome variable

# Training the deep learning model with 5 hidden layers
dlearning.mlp <- h2o.deeplearning(y = y.dep,
                                  x = x.indep,
                                  training_frame = train,
                                  validation_frame = test,
                                  distribution = "bernoulli",
                                  activation = "RectifierWithDropout",
                                  hidden = c(100, 100, 100, 100, 100),
                                  input_dropout_ratio = 0.2,
                                  hidden_dropout_ratios = c(0.5, 0.5, 0.5, 0.5, 0.5),
                                  l2 = 1e-1,
                                  variable_importances = TRUE,
                                  epochs = 100)

# Ranking of predictors in terms of their importance
imp.mlp <- as.data.frame(h2o.varimp(dlearning.mlp))

# Performance on the test set
h2o.performance(dlearning.mlp, valid = TRUE)

# Alternatively, we can use the following code to obtain the probability of M and the predicted class for the test cases
predict.dl <- as.data.frame(h2o.predict(dlearning.mlp, test))
P <- predict.dl$M
C <- predict.dl$predict

mean(C==TE$y)
table(TE$y, C)
auc <- roc(TE$y, P, plot=TRUE)
legend('bottomright', legend=paste('AUC =', round(auc$auc*100, 1)), bty = 'n')
