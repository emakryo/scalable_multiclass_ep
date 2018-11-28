This folder contains the R code of the 4 methods described in the main paper
EP, SEP, VI and GFITC. The code can be evaluated in the MNIST dataset for the 
stochastic version of the methods. For that, the first step is to generate the 
data. For this, you can use the supplied script "download_mnist.R" in the "data" 
folder, which downloads the needed files from http://yann.lecun.com/exdb/mnist/ 
and saves them in a file called "mnist.dat"

To do that, you should go to that folder and run 

R --no-save < download_mnist.R 

That will generate the MNIST data used by the different experiments. 
Batch versions of the methods will run on the Satellite dataset from the UCI 
repository. In the "data" folder we have included a file "satellite.dat", as
well as another file called "folds.dat", that contains 100 splits of the dataset

Then, in the folder "experiments"  there is one sub-folder per method with one
script inside each sub-folder of type "simulate_ ...". Simply run that script 
using 

R --no-save < name_of_the_script

That will simulate the corresponding method and the results will be stored
in the sub-folder results. When launching a batch method the number of the
split needs to be specified, using

R --no-save < name_of_the_script split

Note that several R packages have to be installed for the code to work properly.
You can install them using the typical R package installing tool.

