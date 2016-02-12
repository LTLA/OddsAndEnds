rm logcount.*
bsub -J "counter" -R "rusage[mem=16000]" -n 4 -e "logcount.err" -o "logcount.out" bash count_me.sh
