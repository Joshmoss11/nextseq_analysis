#!/bin/tcsh
#SBATCH --mem=16000
#SBATCH --cpus-per-task=32
#SBATCH --output=run.out
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
if ($#argv != 1) then
        echo "Usage: $0 <Run Name>"
        exit 0
endif
#echo `ls /cs`
#ls /cs/nextseq/
#ls /cs/nextseq/${1}/
#cp /cs/nextseq/${1}/RunInfo.xml ./

#cd /cs/nextseq/${1}/
#ls
#echo Nextseq data can be found in /cs/nextseq/${1} > run.info
#echo converting fastq files
#srun bcl2fastq -i /cs/nextseq/${1}/Data/Intensities/BaseCalls/ -o ./ --create-fastq-for-index-reads --no-lane-splitting --no-bgzf-compression --sample-sheet /cs/icore/joshua.moss/scripts/nextseq/SampleSheet.csv
#rm -rf temp* InterOp Stats RunInfo.xml Reports
#echo done converting fastq files
#echo extracting fastq files
#srun gunzip ./*gz
#echo done extracting fastq files
echo beginning matlab analysis
setenv MATLABPATH '/cs/icore/joshua.moss/scripts/nextseq'
matlab -nodesktop -nosplash -r 'methyl_plasma_pipeline3'
rm -f *fastq *fasta
echo Done analysis!
