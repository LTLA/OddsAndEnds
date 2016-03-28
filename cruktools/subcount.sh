if [[ -e logcount.err || -e logcount.out ]]
then
    rm logcount.*
fi

torun=count_me.sh
if [[ $# -eq 1 ]]
then
    torun=$1
fi

bsub -J "counter" -R "rusage[mem=16000]" -n 4 -e "logcount.err" -o "logcount.out" bash $torun
