set -e
set -u 

# This Bash script is designed to provide support for the ChIA-PET mapping pipeline
# i.e. to go from *.sra files to *.bam files.  It does so in a manner that is
# robust to crashes i.e. can be reloaded at each command if necessary. I've refrained
# from putting in an e.g. Python-mediated pipeline as I'd have to learn all
# the coding conventions, and I don't want to do that until I know which one to
# use (or, have one recommended to me).

###############################################################################
# Making a working space, as well as space to store the finished products
# (i.e. BAM files). We make a very temporary directory for each program to 
# pile working components into. We also need space to store some diagnostic 
# plots.

temp=working
if [ ! -e $temp ]; then
		mkdir $temp
fi
bam=bam
if [ ! -e $bam ]; then
		mkdir $bam
fi
vtmp=localtmp
if [ ! -e $vtmp ]; then
		mkdir $vtmp
fi
pics=pics
if [ ! -e $pics ]; then
		mkdir $pics
fi

# We load some common functions that we need from a BASH mapping library.
# This includes all the logging functions and something to guess the Phred
# offset so that Bowtie2 doesn't throw up.

log=trim.log
source `dirname "${BASH_SOURCE[0]}"`/maplib.sh

###############################################################################

for sra in ${files[@]}; do
		lastcheck=$sra":endsort"
		if [ `check_done $lastcheck $log` -ne 1 ]; then
			continue
		fi
		if [[ $sra =~ \.sra$ ]]; then
			prefix=`basename $sra | sed "s/\.sra$//g"`
		else
			prefix=`basename $sra | sed -E "s/\.(fastq|fq)\.gz$//g"`
		fi

		# Decompressing.
		curjob=$sra":split"
		if [ `check_done $curjob $log` -eq 1 ]; then
			if [[ $sra =~ \.sra$ ]]; then
				fastq-dump -O $temp --split-files $sra 
			else 
				# Just to deal with crazy upload formats.
				alt=`echo $sra | sed -E "s/_1\.((fastq|fq)\.gz)$/_2.\1/"`
				if [[ $alt == $sra ]]; then
					continue
				fi
				zcat $sra > $temp/${prefix}_1.fastq
				zcat $alt > $temp/${prefix}_2.fastq
			fi
			echo $curjob >> $log
		fi

		# Checking the files, if it just got split.
		fastq1=`find $temp | egrep "${prefix}_1\.(fastq|fq)$"`
		fastq2=`find $temp | egrep "${prefix}_2\.(fastq|fq)$"`
		if [[ $fastq1 == "" ]] || [[ $fastq2 == "" ]]; then
			echo "Paired FastQ files not found" >&2
			exit 1
		fi
		prefix1=`echo $fastq1 | sed -r "s/_1\.(fastq|fq)$//"`
		prefix2=`echo $fastq2 | sed -r "s/_2\.(fastq|fq)$//"`
		if [[ $prefix1 != $prefix2 ]]; then
			echo "Prefixes do not match up for the two FastQ files"
			exit 1
		fi

		# Making QC plots (removing the *.zip file that goes with it).
		curjob=$sra":plot"
		if [ `check_done $curjob $log` -eq 1 ]; then
			fastqc $fastq1 -o $pics
			fastqc $fastq2 -o $pics
			rm $pics/*.zip
			echo $curjob >> $log
		fi

		# Linker identification.
		curjob=$sra":delinked"
		tempfix=$temp/$prefix
		if [ `check_done $curjob $log` -eq 1 ]; then
			additional=${additional:=""}
			${HOME}/petchia/scripts/splitLinkers/splitLinkers -a $linkA -b $linkB -1 $fastq1 -2 $fastq2 -o $tempfix $additional
			echo $curjob >> $log
		fi

		# Aligning each linker set with Bowtie2 (removing once done, so it doesn't attempt to realign if one fails.
		curjob=$sra":align"
		if [ `check_done $curjob $log` -eq 1 ]; then
			curphred=`guess_phred $pics $fastq1`
			if [ $curphred -eq 0 ]; then
				curphred=$phred
			fi
			allfiles=(`find $temp -maxdepth 1 -type f | egrep "${tempfix}_[^_]+_[12]\.(fastq|fq)$"`)
			for curf in ${allfiles[@]}; do 
				fout=`echo $curf | sed -r "s/\.(fastq|fq)$/.sam/"`
				bowtie2 -p 8 --reorder --very-sensitive -x $genome -U $curf -S $fout --phred$curphred
				rm $curf
			done
			echo $curjob >> $log
		fi

		# Adding flags and compressing.
		curjob=$sra":recon"
		if [ `check_done $curjob $log` -eq 1 ]; then
			allfiles=(`find $temp -maxdepth 1 -type f | egrep "_[12]\.sam$"`)
			for curf in ${allfiles[@]}; do 
				if [[ $curf =~ _1\.sam$ ]]; then
					additive=64
				elif [[ $curf =~ _2\.sam$ ]]; then
					additive=128
				else
					echo "WHAT?"
					exit 1
				fi
				outcome=`echo $curf | sed "s/\.sam$/.bam/"`
				echo $curf $outcome
				awk -v extra=$additive 'BEGIN { OFS = "\t" } { if (substr($1, 0, 1)!="@") $2=$2+1+extra; print $0 }' $curf | samtools view -bS - > $outcome
				rm $curf
			done
			echo $curjob >> $log
		fi

		# Merging files.
		curjob=$sra":merge"
		if [ `check_done $curjob $log` -eq 1 ]; then
			allfiles=(`find $temp -maxdepth 1 -type f | egrep "_1\.bam$"`)
			for curf in ${allfiles[@]}; do
				matef=`echo $curf | sed "s/_1\.bam$/_2.bam/"`
				finalf=`echo $curf | sed "s/_1\.bam$/.bam/"`
				samtools merge -nf $finalf $curf $matef
				rm $matef $curf
			done
			echo $curjob >> $log
		fi
		
		# Fixing the Mate information (continuing if it's already coordinate sorted).
		curjob=$sra":fixmate"
		if [ `check_done $curjob $log` -eq 1 ]; then
			allfiles=(`find $temp -maxdepth 1 -type f | egrep "\.bam$"`)
			for curf in ${allfiles[@]}; do
				if [ `samtools view -H $curf | grep "SO:coordinate" | wc -l` -ne 0 ]; then
					continue
				fi
				FixMateInformation I=$curf TMP_DIR=$vtmp SO=coordinate VALIDATION_STRINGENCY=SILENT
			done
			echo $curjob >> $log
		fi

		# Marking duplicates (continuing it's already been marked).
		curjob=$sra":dedup"
		if [ `check_done $curjob $log` -eq 1 ]; then
			tempname=`mktemp $temp/blahXXXXX`".bam"
			allfiles=(`find $temp -maxdepth 1 -type f | egrep "\.bam$"`)
			for curf in ${allfiles[@]}; do
				if [ `samtools view -H $curf | egrep "MarkDuplicates" | wc -l` -eq 1 ]; then
					continue
				fi
				MarkDuplicates I=$curf O=$tempname M=$temp/blah.txt TMP_DIR=$vtmp AS=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
				mv $tempname $curf
			done
			echo $curjob >> $log
		fi

		# Resorting by name.
		allfiles=(`find $temp -maxdepth 1 -type f | egrep "\.bam$"`)
		for curf in ${allfiles[@]}; do
			leftover=`basename $curf | sed "s/\.bam$//"`
			samtools sort -n $curf $bam/$leftover
			rm $curf
		done
		echo $lastcheck >> $log

		# Cleaning out.
		rm $temp/*
done

###############################################################################
# Erasing the leftover directories, once it's clear that everything has been
# processed successfully (no need to hold onto them).

rmdir $temp
rm $vtmp/*
rmdir $vtmp
rm $log   

###############################################################################
###############################################################################
###############################################################################
# Finish.


