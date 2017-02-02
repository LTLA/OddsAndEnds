set -e 
set -u

# Given an input CRAM file, this script will generate one or two FASTQ files.
# If you give it two files, it'll assume that the reads are paired and behave appropriately.

if [ $# -eq 2 ]
then
    paired=0
elif [ $# -eq 3 ]
then
    paired=1
else
    echo "$0 CRAM FASTQ [MATE]"
    exit 1
fi

prefix=$(basename $1 | sed -r "s/\\.cram$//")
odir=$(dirname $2)
workingfix=${odir}/tempcram_${prefix}
working=${workingfix}.bam

echo samtools collate $1 ${workingfix}
if [ $paired -eq 0 ]
then
    echo samtools fastq -F 2304 ${working} -s ${2}
else
    echo samtools fastq -F 2304 ${working} -1 ${2} -2 ${3}
fi

rm ${working}
exit 0
