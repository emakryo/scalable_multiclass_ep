##########################################################
# This script downloads and prepares the needed files 
#           to work with the MNIST dataset
##########################################################

library("R.utils")

if (! file.exists("train-images-idx3-ubyte")) {
    if (! file.exists("train-images-idx3-ubyte.gz"))
        download.file("http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz", "train-images-idx3-ubyte.gz")
    gunzip("train-images-idx3-ubyte.gz")
}

if (! file.exists("t10k-images-idx3-ubyte")) {
    if (! file.exists("t10k-images-idx3-ubyte.gz"))
        download.file("http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz", "t10k-images-idx3-ubyte.gz")
    gunzip("t10k-images-idx3-ubyte.gz")
}

if (! file.exists("train-labels-idx1-ubyte")) {
    if (! file.exists("train-labels-idx1-ubyte.gz"))
        download.file("http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz", "train-labels-idx1-ubyte.gz")
    gunzip("train-labels-idx1-ubyte.gz")
}
   
if (! file.exists("t10k-labels-idx1-ubyte")) {
    if (! file.exists("t10k-labels-idx1-ubyte.gz")) 
        download.file("http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz", "t10k-labels-idx1-ubyte.gz")
    gunzip("t10k-labels-idx1-ubyte.gz")
}


to.read <- file("t10k-images-idx3-ubyte", "rb")
readBin(to.read, integer(), n=4, endian="big")

test_data <- matrix(0, 1e4, 28 * 28)

for (i in 1 : 1e4) {
    
    if (i %% 100 == 0) cat(".")
    
    m <- matrix(as.integer(readBin(to.read, raw(), size = 1, n = 28 * 28, endian="big")), 28, 28)
    test_data[ i, ] <- c(m)
}

cat("\n")

to.read <- file("t10k-labels-idx1-ubyte", "rb")
readBin(to.read, integer(), n=2, endian="big")

test_labels <- as.integer(readBin(to.read, raw(), size = 1, n = 1e4, endian="big"))

###

to.read <- file("train-images-idx3-ubyte", "rb")
readBin(to.read, integer(), n=4, endian="big")

train_data <- matrix(0, 6e4, 28 * 28)

for (i in 1 : 6e4) {
    
    if (i %% 100 == 0) cat(".")
    
    m <- matrix(as.integer(readBin(to.read, raw(), size = 1, n = 28 * 28, endian="big")), 28, 28)
    train_data[ i, ] <- c(m)
}

cat("\n")

to.read <- file("train-labels-idx1-ubyte", "rb")
readBin(to.read, integer(), n=2, endian="big")

train_labels <- as.integer(readBin(to.read, raw(), size = 1, n = 6e4, endian="big"))

data <- list(Xtr = train_data, Ytr = train_labels, Yts = test_labels, Xts = test_data)

save(data, file = "mnist.dat")


