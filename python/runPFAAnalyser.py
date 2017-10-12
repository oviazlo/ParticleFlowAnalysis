#!/bin/python
import subprocess

fileDir = '/ssd/viazlo/data/CLIC_o3_v13_ILCSoft-2017-08-23_gcc62_photons_v10_files/'
filePrefix = 'CLIC_o3_v13_ILCSoft-2017-08-23_gcc62_photons_v10'

energyToUse = [10]
thetaToUse = [30,31,32,33,34,35,36,37,38,39,40]
bashCommand = 'bash /afs/cern.ch/user/v/viazlo/SCRIPTS/executeParallelThreads.sh'

counter = 0

for iE in energyToUse:
	for iTheta in thetaToUse:
		bashCommand = bashCommand + ' "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/data/PhotonEfficiency ' + fileDir + filePrefix + '_E' + str(iE) + '_theta' + str(iTheta) + '*"'
		counter = counter + 1
		if (counter > 7):
			#  process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
			#  output, error = process.communicate()
			print (bashCommand)
			bashCommand = 'bash ~/SCRIPTS/executeParallelThreads.sh'
			counter = 0
if (counter != 0):
	#  process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        #  output, error = process.communicate()
	print (bashCommand)
