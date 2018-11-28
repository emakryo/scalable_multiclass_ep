source("../../../Rcode/ep_code/epMGPC_stochastic.R")

set.seed(0)

load("../../../data/mnist.dat")

fileName <- "./results/mnist-01-"

s <- sample(1 : nrow(data$Xtr))[ 1 : 6e4 ]

data$Xtr <- data$Xtr[ s , ]
data$Ytr <- as.factor(data$Ytr[ s ])
data$Yts <- as.factor(data$Yts)

# Data split
X_train <- data$Xtr
Y_train <- as.factor(data$Ytr)

X_test <- data$Xts
Y_test <- as.factor(data$Yts)

m <- 200

# Normalize

X_train <- (X_train / 255)

X_test <- (X_test / 255)

D <- as.matrix(dist(X_train[ sample(1 : nrow(X_train), 1e3), ]))
l <- median(D[ upper.tri(D) ])
log_l = matrix(log(l), 10, ncol(X_train))

REPORT <- TRUE
CONT <- 1

t0 <- proc.time()

ret <- epMGPCInternal(X_train, Y_train, m, m, X_test = X_test, Y_test = Y_test, optim = "adam", print_interval = 10, log_l = log_l)

t_after <- proc.time()

