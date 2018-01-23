
'''

Python script to submit CRAB3 jobs processing several Monte Carlo samples

Run the following commands from the parent directory: 
(you might want to double-check paths and copy scripts to the correct locations)

1. Submit CRAB3 jobs
    cmsenv && source /cvmfs/cms.cern.ch/crab3/crab.sh && python MCTuples/crabConfig_mc.py
If the jobs crash (error 8002), try checking that in the CMSSW config 'runOnVM' is False!

2. Merge the tuples of each sample: get one file per pthat interval
    bsub -q 1nh -J job1 < MCTuples/mergeMCTuples.sh
3. Create the final Monte Carlo tuple file
(the ROOT script does the merging, while also computing: weight = cross-section / num of events):
    bsub -q 1nh -J job1 < MCTuples/createMCTuple.sh 

'''


from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

# Global CRAB3 settings
config.General.requestName = 'OpenDataTree_QCDFlat'
config.General.workArea = 'crab_projects_retry'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'OpenDataTreeProducerOptimized_mcPAT_2012_cfg.py'
config.JobType.outputFiles = ['ak5PFJetsParams0.npy','ak7PFJetsParams0.npy']
#config.JobType.inputFiles = ['/uscms_data/d2/aperloff/YOURWORKINGAREA/MLJEC/CMSSW_5_3_32/src/fastjet-contrib/lib/']
config.JobType.sendExternalFolder = True

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' #'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/lpcjme/noreplica/'
config.Data.publication = False
config.Data.ignoreLocality = True

config.Data.totalUnits = -1         # appr. number of events
config.Data.unitsPerJob = 2     # appr. events per job

config.Site.storageSite = 'T3_US_FNALLPC'


# The following chunk was taken from:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRABClientLibraryAPI?rev=29
if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    # Helper functions
    def submit(config):
        try:
            #crabCommand('submit', config = config)
            crabCommand('submit', config = config, dryrun = True) # Debug purposes
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    def extractName(dataset):
        # Return shorter dataset name, e.g. "pthat_15to30"
        import re
        m = re.search('QCD_Pt-(\d+to\d+)_TuneZ2star_Flat_8TeV_pythia6', dataset.split('/')[1])
        return "pthat_" + m.group(1)


    # Loop over MC datasets
    for dataset in ['/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-NoPileUp_START53_V7N-v1/AODSIM']:

        #config.Data.inputDataset = dataset
        config.Data.userInputFiles = open('CMS_MonteCarlo2012_Summer12_DR53X_QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_AODSIM_NoPileUp_START53_V7N-v1_combined_file_index.txt').readlines()
        config.Data.outputPrimaryDataset = 'QCD_Pt-15to3000_Flat_8TeV_NoPileUp_v1'
        config.General.requestName = extractName(dataset)


        # We use a multiprocessing trick explained here:
        # https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ?rev=56#Multiple_submission_fails_with_a
        # Reason: calling repeatedly 'submit(config)' gives an error (possibly related to the PAT modules)
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

