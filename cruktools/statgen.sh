# This script compiles mapping and quality statistics for all BAM files in a directory.

set -e 
set -u

tmpdir=$(mktemp -d --tmpdir=.)
temp=$tmpdir/stat.log
output=all_qual.tsv
echo -e "Sample\tTotal\tMapped\tPaired\tProperPair\tBothMapped\tInterChr\tFiltered\tDedupped\tBasicStat\tPerBaseQual\tPerTileQual\tPerSeqQual\tPerBaseSeq\tPerSeqGC\tPerBaseN\tSeqLen\tSeqDup\tOverSeq\tAdaptorContent\tKmerContent" > $output

for x in $(ls bam | grep "\\.bam$")
do
	curfile=bam/$x
	prefix=$(echo $x | sed "s/\\.bam$//")
	
	# Deparsing samtools flagstat.
	samtools flagstat $curfile > $temp
	total=$(cat $temp | grep "in total (" | cut -f1 -d " ")
	mapped=$(cat $temp | grep "mapped (" | cut -f1 -d " ")
	paired=$(cat $temp | grep "paired in sequencing" | cut -f1 -d " ")
	proper=$(cat $temp | grep "properly paired" | cut -f1 -d " ")
	bothmap=$(cat $temp | grep "with itself and mate mapped" | cut -f1 -d " ")
	oddmap=$(cat $temp | grep "with mate mapped to a different chr$" | cut -f1 -d " ")

	# Probably could cut half these steps with flagstat, but more effort to parse.
	filtered=$(samtools view -c -q 10 -F 4 $curfile)
	dedupped=$(samtools view -c -q 10 -F 1028 $curfile)

    # Unpacking quality score files (possibly many if paired-end).
    echo $prefix
    for y in $(ls qual | egrep "${prefix}(_[12])?.zip")
    do
    	unzip -qq qual/$y -d $tmpdir

        # Parsing summaries
        sumfile=$tmpdir/*/summary.txt
        basicstats=$(cat $sumfile | grep "Basic Statistics" | cut -f1)
        perbasequal=$(cat $sumfile | grep "Per base sequence quality" | cut -f1)
        pertilequal=$(cat $sumfile | grep "Per tile sequence quality" | cut -f1)
        perseqqual=$(cat $sumfile | grep "Per sequence quality scores" | cut -f1)
        perbaseseq=$(cat $sumfile | grep "Per base sequence content" | cut -f1)
        perseqgc=$(cat $sumfile | grep "Per sequence GC content" | cut -f1)
        perbaseN=$(cat $sumfile | grep "Per base N content" | cut -f1)
        seqlen=$(cat $sumfile | grep "Sequence Length Distribution" | cut -f1)
        seqdup=$(cat $sumfile | grep "Sequence Duplication Levels" | cut -f1)
        overseq=$(cat $sumfile | grep "Overrepresented sequences" | cut -f1)
        adapt=$(cat $sumfile | grep "Adapter Content" | cut -f1)
        kmer=$(cat $sumfile | grep "Kmer Content" | cut -f1)

        rm -r $tmpdir/*
    	echo -e "$prefix\t$total\t$mapped\t$paired\t$proper\t$bothmap\t$oddmap\t$filtered\t$dedupped\t$basicstats\t$perbasequal\t$pertilequal\t$perseqqual\t$perbaseseq\t$perseqgc\t$perbaseN\t$seqlen\t$seqdup\t$overseq\t$adapt\t$kmer" >> $output
        
        total='-'
        mapped='-'
        paired='-'
        proper='-'
        bothmap='-'
        oddmap='-'
        filtered='-'
        dedupped='-'
    done
done
rm -r $tmpdir
