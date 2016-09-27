#!/bin/tcsh
#SBATCH --job-name=nextseq
#SBATCH --output=run.out
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --nodes=2
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
if ($#argv != 1) then
        echo "Usage: $0 <Run Name>"
        exit 0
endif
#echo `ls /cs`
#ls /cs/nextseq/
#ls /cs/nextseq/${1}/
cp /cs/nextseq/${1}/RunInfo.xml ./

#cd /cs/nextseq/${1}/
#ls
echo Nextseq data can be found in /cs/nextseq/${1} > run.info
echo converting fastq files
bcl2fastq -i /cs/nextseq/${1}/Data/Intensities/BaseCalls/ -o ./ --create-fastq-for-index-reads --no-lane-splitting --no-bgzf-compression --sample-sheet /cs/icore/joshua.moss/scripts/nextseq/SampleSheet.csv
rm -rf temp* InterOp Stats RunInfo.xml Reports
echo done converting fastq files
echo extracting fastq files
gunzip ./*gz
echo done extracting fastq files
echo beginning matlab analysis
mkdir -p /tmp/$USER/$SLURM_JOB_ID
setenv MATLABPATH '/cs/icore/joshua.moss/scripts/nextseq'
matlab -nodisplay -r  'methyl_plasma_pipeline4'
rm -rf *fastq *fasta /tmp/$USER/$SLURM_JOB_ID
echo Done analysis!
