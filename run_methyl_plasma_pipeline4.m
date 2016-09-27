disp('JOB ID IS')
disp(getenv('SLURM_JOB_ID'))
%oldF = cd(strcat('/cs/icore/',getenv('USER')))
pc = parcluster('local');
pc.NumWorkers=64;
pc.JobStorageLocation = strcat('/tmp/',getenv('USER'),'/', getenv('SLURM_JOB_ID'));

parpool(pc, pc.NumWorkers);
%cd(oldF)
methyl_plasma_pipeline4
