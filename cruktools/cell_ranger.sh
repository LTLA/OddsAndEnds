set -e
set -u

if [ $# -le 2 ]
then
    echo "$0 <barcode> <annotation> [OPTIONS]"
    exit 1
fi
bc=$1
anno=$2
shift
shift
otheropts="$@"

# Setting up the working directory.
workdir=working_${bc}
if [ ! -e ${workdir} ]
then
    mkdir ${workdir}
fi
StripeThis ${workdir}

# Setting up the final output directory.
if [ ! -e results ]
then
    mkdir results
fi

# Unpacking the files into 'working' and setting up the output directory.
if [ ! -e ${bc} ]
then 
    for f in $(ls fastq | grep "${bc}.*tar$")
    do
        tar xf fastq/${f} -C ${workdir}
    done
    mkdir ${bc}
    StripeThis ${bc}
fi

# Setting up the output directory.
newdir=results/${bc}
if [ ! -e ${newdir} ]
then
    mkdir ${newdir}
fi

# Executing cellranger
/home/mib-cri/software/10Xgenomics/cellranger-1.2.0/cellranger count \
    --id="${bc}" \
    --transcriptome="${anno}" \
    --fastqs="${workdir}" \
    --sample="${bc}" \
    ${otheropts}

# Pulling out the analysis files.
cp -L ${bc}/outs/filtered_gene_bc_matrices/mm10/* ${newdir}
cp -L ${bc}/outs/web_summary.html ${newdir}
cp -L ${bc}/outs/metrics_summary.csv ${newdir}

# Cleaning out the files.
rm -r ${workdir}
rm -r ${bc}

