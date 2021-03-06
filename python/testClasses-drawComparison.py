from __future__ import print_function
import unittest
from classesForDrawComparison import helperClass, nLegendCaptions, notUsedAttribute, noMandatoryAttribute, pythonEvalFuncError
import sys
import yaml # WARNING use this environment ~/env/setupPyTools27.env to enable yaml package!!!
import ROOT

# Main body of the test suite starts here
class QueueTest(unittest.TestCase):

    def readConfig(self):
        with open('/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/FOR_TEST_DO_NOT_REMOVE.yml', 'r') as yamlFile:
            globalCfg = yaml.load(yamlFile)
            #  print ("[info]\t Read yaml file: \n\t\t%s" % (yamlFile))

            if (globalCfg.get("default") is not None):
                defaultCfg = globalCfg.get("default")
                for cfgIterator in globalCfg:
                    cfg = globalCfg[cfgIterator]
                    if (cfg == defaultCfg):
                        continue
                    for cfgIt in defaultCfg:
                        if (cfg.get(cfgIt) is None):
                            cfg[cfgIt] = defaultCfg.get(cfgIt)
                globalCfg.pop("default")
            self.globalCfg = globalCfg

    def setUp(self):
        self.q = helperClass()
        self.readConfig()

    def testConstructor(self):
        self.assertTrue(self.globalCfg.get('default') is None)
        self.assertTrue(len(self.globalCfg) == 9)

    def testGetProcessedHists(self):
        cfg = self.globalCfg['thetaRes_DR11_FCCee_vs_CLIC_vsCosTheta']
        hists = self.q.getProcessedHists(cfg)
        self.assertEqual(len(hists), 4)
        self.assertEqual(list(x.GetMarkerStyle() for x in hists), [4,5,6,ROOT.kOpenCircle]  )
        self.assertEqual(list(x.GetMarkerColor() for x in hists), [777,888,999,ROOT.kBlack]  )

    def testGetLegend(self):
        cfg = self.globalCfg['thetaRes_DR11_FCCee_vs_CLIC_vsCosTheta']
        #  print(cfg)
        hists = self.q.getProcessedHists(cfg)
        tmpCfg = dict(cfg)
        self.assertRaises(nLegendCaptions, self.q.getLegend,hists,tmpCfg)
        cfg['legTitle'] += ['#sqrt{s} = 365 GeV, DR09']
        leg = self.q.getLegend(hists,cfg)
        # second call of the functino getLegend:
        self.assertTrue(self.q.getLegend(hists,cfg))
        #  print(cfg)
        
    def testGetTextLabels(self):
        cfg = dict(self.globalCfg['testTextLabels'])
        cfg.pop('histName')
        self.assertRaises(noMandatoryAttribute,self.q.getProcessedHists,cfg)

        cfg = dict(self.globalCfg['testTextLabels2'])
        hists = self.q.getProcessedHists(cfg)
        self.assertRaises(notUsedAttribute,self.q.getTextLabels,hists,cfg)

        # attempts to use eval() with wrong input type
        cfg = dict(self.globalCfg['testLatexLabels1'])
        hists = self.q.getProcessedHists(cfg)
        self.assertRaises(pythonEvalFuncError,self.q.getLatexLabels,hists,cfg)

        # attempts to use eval() for set-function with multiple arguments (not implemented yet)
        cfg = dict(self.globalCfg['testLatexLabels2'])
        hists = self.q.getProcessedHists(cfg)
        self.assertRaises(pythonEvalFuncError,self.q.getLatexLabels,hists,cfg)

    def testTF1_nFuncs(self):
        cfg = dict(self.globalCfg['testTF1_1'])
        funcs = self.q.getTF1s(cfg)
        self.assertTrue(len(funcs)==2)

    def testTF1_wrongFuncName(self):
        cfg = dict(self.globalCfg['testTF1_2'])
        funcs = self.q.getTF1s(cfg)
        self.assertTrue(len(funcs)==1)

# ----------------------------------------------------------------------

# Unittest does not expect to find the extra parameters we pass (time
# limit, implementation). Remove them before running unittest.
sys.argv = sys.argv[:1]

try:
    # unittest shuts down the interpreter when it finishes the
    # tests. We want to delay the exit, in order to display some
    # timing information, so we handle the SystemExit exception.
    unittest.main()
except SystemExit:
    pass

