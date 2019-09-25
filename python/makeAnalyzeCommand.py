#!/bin/python

#  availableParticleTypes = ["e-","mu-","gamma","pi-"]
availableParticleTypes = ["SinglePart_gamm"]
availableParticleEnergies = [5,10,20,50,100]

scriptPath = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/build_20jun2019/PhotonECAL/StudyElectronPerformance"
detModel = "FCCee_o1_v03"
softVer = "ILCSoft-2018-04-26_gcc62"
dataPath = "/eos/experiment/clicdp/grid/ilc/user/o/oviazlo/"
productionTag = "ConformalTr_run1_50mm"


for particleType in availableParticleTypes:
    outDir = particleType
    f = open(outDir + '/runIt.sh', 'w')
    for iE in availableParticleEnergies:
      runCommand = 'nohup ' + scriptPath + ' -f "' + dataPath + detModel + '_' + softVer + '_' + particleType + '_E' + str(iE) + '_' + productionTag + '_files/' + detModel + '_' + softVer + '_' + particleType + '_E' + str(iE) + '_' + productionTag  + '_*" --energy ' + str(iE-1) + ' ' + str(iE+1) + ' --theta 50 130 &'
      print (runCommand)
      f.write(runCommand+'\n')
    f.close()



