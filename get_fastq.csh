f ($#argv != 1) then
        echo "Usage: $0 <Run Name>"
        exit 0
endif
cp /cs/nextseq/${1}/RunInfo.xml ./
bcl2fastq -i /cs/nextseq/${1}/Data/Intensities/BaseCalls/ -o fastq_files --create-fastq-for-index-reads --no-lane-splitting

