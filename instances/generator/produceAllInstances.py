# this script has to be execute with this command line:
# python produceAllInstances.py

from commands import getstatusoutput
import os.path

# the directory to put the instances
dirInstances = ".."

# the directory where the generator can be found
dirGenerator = "."

# list of instance parameters to produce
# of course you can change this list
listRho = [-0.9, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 0.9]
listM   = [2, 3, 5]
listN   = [18, 32, 64, 128]
listK   = [2, 4, 6, 8, 10]

# number of instance per tuple of parameters
nbInst = 30

# create the directory of instances if necessary
if not os.path.isdir(dirInstances):
    os.mkdir(dirInstances)

# Go to produce all the instances. 
# It can take a lot of time (hours) if there is a lot of large instances !
for M in listM:
    for N in listN:
        for K in listK:
            if K < N:
                for r in listRho:
                    for n in range(nbInst):
                        name = "rmnk_" + str(r) + "_" + str(M) + "_" + str(N) + "_" + str(K) + "_" + str(n) + ".dat"
                        
                        if r > (-1.0 / (M - 1)):
                            command_line = "R --slave --no-restore --file=" + dirGenerator + "/rmnkGenerator.R --args " + str(r) + " " + str(M) + " " + str(N) + " " + str(K) + " " + str(n) + " " + dirInstances + "/" + name
                            # you can use this command if you need:
                            # command_line = dirGenerator + "/rmnkGenerator.R " + str(r) + " " + str(M) + " " + str(N) + " " + str(K) + " " + str(n) + " " + dirInstances + "/" + name
                            res = getstatusoutput(command_line)
