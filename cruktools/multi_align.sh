set -e
set -u

runfailed=${runfailed:=0}
extra=${extra:=""}
if [ -z ${ispet+x} ]
then
    echo "Error: must set PET specification" >&2
    exit 1
fi
if [ -z ${genome+x} ]
then
    echo "Error: must set genome index" >&2
    exit 1
fi
if [ -z ${fastq+x} ]
then
    echo "Error: must set FASTQ array" >&2
    exit 1
fi

# Setting up folders to put stuff in (do it here, otherwise race conditions cause jobs to fail when object is constructed in another thread just after the test).
if [ ! -e qual ]
then
    mkdir qual
fi

if [ ! -e bam ]
then
    mkdir bam
    lfs setstripe -s 64m -c -1 bam/ # Adding stripe settings.
fi

if [ ! -e logs ]
then
    mkdir logs
fi

jname=Align$RANDOM
location=$( dirname $BASH_SOURCE )
for x in ${fastq[@]}
do
    subsec=$(basename $x)

    # Processing for FASTQ files
    if [[ $subsec =~ "(fastq|fq)" ]]
    then
        subsec=$(echo $subsec | sed -r "s/\\.(fastq|fq)(\\.gz)?$//")
        supercmd="eval ${location}/solo_align.sh -f $x -i $genome ${extra}" # Added eval in case 'extra' has quotes.
        if [[ $ispet -ne 0 ]]
        then
            if [[ $x =~ "1\\.(fastq|fq)" ]] 
            # Skipping if it's not the first read.
            then
                subsec=$(echo $subsec | sed -r "s/_?(p|R)?1$//")
                mate=$(echo $x | sed -r "s/1\\.(fastq|fq)/2.\1/")
                supercmd="${supercmd} -m ${mate}"
            else
                continue
            fi
        fi
        supercmd="${supercmd} -p ${subsec}"

    # Processing for CRAM files; first to BAM, then to (paired-end) FASTQ
    elif [[ $subsec =~ "cram" ]]
    then 
        subsec=$(echo $subsec | sed -r "s/\\.cram$//")
        workingfix=bam/tempcram_${subsec}
        working=${workingfix}.bam
        aligncmd="eval ${location}/solo_align.sh -i ${genome} -p ${subsec} ${extra}"
        if [[ $ispet -eq 0 ]]
        then
            ref=bam/temp_${subsec}.fastq
            supercmd="set -e; set -u; bash ${location}/cram2fastq.sh $x ${ref}; ${aligncmd} -f ${ref}; rm ${ref}"
        else 
            first=bam/temp_${subsec}_1.fastq
            mate=bam/temp_${subsec}_2.fastq
            supercmd="set -e; set -u; bash ${location}/cram2fastq.sh $x ${first} ${mate}; ${aligncmd} -f ${first} -m ${mate}; rm ${first} ${mate}"
        fi
    fi
 
    # Only running failed jobs, if requested.
    if [ -e logs/${subsec}.log ] && [ $runfailed -eq 1 ]
    then
        continue
    fi        

    # Deleting existing logs
    for f in $(ls logs | grep "^${subsec}\\.")
    do
        rm logs/$f
    done

    # Adding a job and polling. 
    bsub -J "$jname" -R "rusage[mem=16000]" -n 1 -e "logs/${subsec}.err" -o "logs/${subsec}.out" ${supercmd}
    while [ $( bjobs -J $jname | wc -l ) -gt 11 ]; do sleep 10; done
done
    
# If the 'log' has been reported, we delete the 'err' and 'out' files for that run.
while [ $( bjobs -J $jname | wc -l ) -gt 0 ]; do sleep 10; done
for x in $( ls logs | grep "\\.log$" )
do
    prefix=$( echo $x | sed "s/\\.log$//" )
    if [ -e logs/${prefix}.err ]
    then
        rm logs/${prefix}.err
        rm logs/${prefix}.out
    fi
done

