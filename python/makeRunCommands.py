#!/bin/python
import sys
energy = [1,2,5,10,20,50,100]
particleType = "pi-"
dirPrefix = "/ssd/viazlo/data/FCCee_o5_v04_ILCSoft-2017-07-27_gcc62_"
dirPostfix = "IsoGun_all_files"
filePrefix = "FCCee_o5_v04_ILCSoft-2017-07-27_gcc62_"

if __name__ == "__main__":

	if (len(sys.argv)>=2):
		particleType = sys.argv[1]

	for iE in energy:
		print ( '/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/build/PhotonECAL/AnalyzeIsotropPhotonSample -f "%s%s%s/%s*E%i_*" --energy %f %f ' % (dirPrefix,particleType,dirPostfix,filePrefix,iE,0.99*iE,1.01*iE ) )
