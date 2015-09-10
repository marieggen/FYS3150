import subprocess
import os
from datetime import datetime

def runInTerminal(cmd):
	print "Running command: ", cmd
	subprocess.call(cmd, shell=True)

compileCommand = "g++ -o project_1_functions.x project_1_functions.cpp"
programName = "./project_1_functions.x"
dataFilePrefix = "data_"
errorFilePrefix = "error_"

runInTerminal(compileCommand)

for n in [1e1, 1e2, 1e3, 1e4, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 1e6]:
	dataFileName = "%s%d.txt" % (dataFilePrefix, n)
	errorFileName = "%s%d.txt" % (errorFilePrefix, n)
	runInTerminal("%s %s %s %d" % (programName, dataFileName, errorFileName, n))