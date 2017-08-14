#!/bin/bash
module load qbic/anaconda2/2.1.0
module load qbic/openms/2.2-14bb753

workflowDir=$(cat wfdir)
#parse using CTDopts and run workflow
python runWorkflow.py $workflowDir
cp wfdir wfdir2