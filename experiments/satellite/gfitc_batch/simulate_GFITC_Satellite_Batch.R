source("../../../Rcode/gfitc_code/epMGPC.R")

load("../../../data/satellite.dat")
load("../../../data/folds.dat")


fileName <- "./results/"


# Fold to train/evaluate

#j <- 1
j <- as.integer(commandArgs()[3])

set.seed(j)


# Train the model for each fold
#   1. Classification error in test 
#   2. Test log likelihood
  
n <- ncol(folds$itrain)
m <- 100
  
X_train <- data$X[ folds$itrain[ j, ] , ]
Y_train <- data$Y[ folds$itrain[ j, ] ]

means <- apply(X_train, 2, mean)
sds <- apply(X_train, 2, sd)
sds[ sds == 0 ] <- 1.0

X_train <- (X_train - matrix(means, nrow(X_train), ncol(X_train), byrow = TRUE)) / 
  matrix(sds, nrow(X_train), ncol(X_train), byrow = TRUE)

# Heuristic for the length scale

M <- as.matrix(dist(X_train))
l <- median(M[upper.tri(M)])

X_test <- data$X[ folds$itest[ j, ] , ]
Y_test <- data$Y[ folds$itest[ j, ] ]

X_test <- (X_test - matrix(means, nrow(X_test), ncol(X_test), byrow = TRUE)) / 
  matrix(sds, nrow(X_test), ncol(X_test), byrow = TRUE)

Xbar <- list()
for (i in 1 : length(levels(Y_train)))
	Xbar[[ i ]] <- X_train[ sample(1 : nrow(X_train),  m), ]

REPORT <- FALSE

t0 <- proc.time()

ret <- epMGPCExternal_lBFGS(X_train, Y_train, m, Xbar_ini = Xbar, l = matrix(l, length(levels(Y_train)), ncol(X_train)))

t_after <- proc.time()

pred_train <- predictMGPC(ret, X_train)
pred_test <- predictMGPC(ret, X_test)

err_train <- mean(pred_train$labels != Y_train) 
neg_log_l_train <- -mean(log(pred_train$prob[ cbind(1:nrow(X_train), as.integer(Y_train))]))
err_test <- mean(pred_test$labels != Y_test) 
neg_log_l_test <- -mean(log(pred_test$prob[ cbind(1:nrow(X_test), as.integer(Y_test))]))

write.table(t(c(t_after - t0)), file = paste(fileName,"time_", j, ".txt", sep=""), col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(err_train, file = paste(fileName,"train-err_", j, ".txt", sep=""), col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(neg_log_l_train, file = paste(fileName,"train-ll_", j, ".txt", sep=""), col.names = FALSE, row.names = FALSE, append = TRUE) 
write.table(err_test, file = paste(fileName,"test-err_", j, ".txt", sep=""), col.names = FALSE, row.names = FALSE, append = TRUE) 
write.table(neg_log_l_test, file = paste(fileName,"test-ll_", j, ".txt", sep=""), col.names = FALSE, row.names = FALSE, append = TRUE) 


