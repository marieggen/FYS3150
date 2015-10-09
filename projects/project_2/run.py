import subprocess
import os
from datetime import datetime

def runInTerminal(cmd):
	print "Running command: ", cmd
	subprocess.call(cmd, shell=True)

# compileCommand = "g++ -std=c++11 -o project_1_functions.x project_1_functions.cpp"
compileCommand = "g++ -o jacobi_tridiag.x jacobi_tridiag.cpp -larmadillo -llapack -lblas"
programName = "./jacobi_tridiag.x"
eigvalues = "eigval_"
eigvectors = "eigvec_"
n = 300

runInTerminal(compileCommand)

o = 20
#0.01, 0.5, 1.0, 5.0
for omega_r in [0.05]:
	eigval = "%s%d%d.txt" % (eigvalues, n, o)
	eigvec = "%s%d%d.txt" % (eigvectors, n, o)
	runInTerminal("%s %s %s %d %f" % (programName, eigval, eigvec, n, omega_r))
	o = o+1
