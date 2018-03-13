#!/bin/python

#  availableParticleTypes = ["e-","mu-","gamma","pi-"]
availableParticleTypes = ["e-"]
availableParticleEnergies = [1,2,5,10,20,50,100]

scriptPath = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/build/PhotonECAL/AnalyzeIsotropPhotonSample"
detModel = "FCCee_o1_v02"
softVer = "ILCSoft-2017-12-21_gcc62"
dataPath = "/ssd/viazlo/data/"
productionTag = "IsoGun_TruthTracking_prod3"


for particleType in availableParticleTypes:
    outDir = particleType
    f = open(outDir + '/runIt.sh', 'w')
    for iE in availableParticleEnergies:
        runCommand = scriptPath + ' -f "' + dataPath + detModel + '_' + softVer + '_' + particleType + productionTag + '_files/' + detModel + '_' + softVer + '_' + particleType + productionTag + '_E' + str(iE) + '_*" --energy ' + str(iE-1) + ' ' + str(iE+1) + ' --theta 9 171'
        print (runCommand)
        f.write(runCommand+'\n')
    f.close()



