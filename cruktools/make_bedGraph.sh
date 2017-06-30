set -e 
set -u

####################################################
# Setting up the options.

discard=0
retain=0
minqual=0
strsplit=0
header=""

while getopts "F:f:q:sh:" opt; do
    case $opt in
        F) # Reads to discard
            discard=$OPTARG
            ;;
        f) # Do genomic alignment.
            retain=$OPTARG
            ;;
        q) # Minimum quality
            minqual=$OPTARG
            ;;
        s) # Split by strand 
            strsplit=1
            ;;
        h) # bedgraph header details
            header=$OPTARG
            ;;
    esac
done

samtools_opts=""
if [[ $minqual -gt 0 ]]
then
    samtools_opts="${samtools_opts} -q $minqual"
fi
if [[ $discard -ne 0 ]]
then
    samtools_opts="${samtools_opts} -F $discard"
fi
if [[ $retain -ne 0 ]]
then
    samtools_opts="${samtools_opts} -f $retain"
fi

shift "$((OPTIND - 1))"
curfile=$1
oprefix=$2

####################################################
# Setting up the intermediate values. 

prefix=$(basename ${curfile} | sed "s/.bam$//")

genfile=${prefix}_genome.txt
samtools view -H $curfile | grep "@SQ" | cut -f2,3 | sed -r "s/[A-Z]+://g" > ${genfile}

bedfile=${prefix}_temp.bed
echo ${samtools_opts}
samtools view $curfile ${samtools_opts} -b | bedtools bamtobed > ${bedfile}

depth=$(wc -l ${bedfile} | cut -f1 -d " ")

####################################################
# Converting the file.

if [[ $strsplit -eq 0 ]]
then
    newname=${oprefix}.bedgraph
    echo "track type=bedGraph name=\"${oprefix}\" ${header}" > $newname
    bedtools genomecov -i ${bedfile} -bg -g ${genfile} -scale $(echo "scale=5;1000000.0/$depth" | bc) >> $newname
else
    for mode in P N 
    do
        oprefix_str=${prefix}_${mode}
        bedfile_str=${oprefix_str}_temp.bed
        if [[ $mode == "P" ]]
        then
            cat ${bedfile} | grep "+" > ${bedfile_str}
        else
            cat ${bedfile} | grep "-" > ${bedfile_str}
        fi

        newname_str=${oprefix_str}.bedgraph
        echo "track type=bedGraph name=\"${oprefix_str}\" ${header}" > $newname_str
        bedtools genomecov -i ${bedfile_str} -bg -g ${genfile} -scale $(echo "scale=5;1000000.0/$depth" | bc) >> $newname_str
    done

    rm ${bedfile_str}
fi

####################################################
# Cleaning out the temporaries.

rm ${genfile} ${bedfile}


