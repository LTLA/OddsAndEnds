set -u
set -e

# This Bash script is designed to provide support for general mapping i.e. to
# go from *.sra/*.fastq.gz files to *.bam files.  It does so in a manner that
# is robust to crashes i.e. can be reloaded at each command if necessary. It 
# will also sort and dedup them (for PE, they can be resorted by name after
# dedupping).

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
# This includes all the logging functions and something to guess the Phred offset.

log=trim.log
source `dirname "${BASH_SOURCE[0]}"`/maplib.sh

###############################################################################
# Now, proceeding to iterate. We make sure we pull out both of the resulting
# *.fastq files from the single *.sra file by using some trickery (i.e. adding
# ...1.fastq and ...2.fastq to the end of each).

for curfile in ${files[@]}; do
	lastcheck=$curfile":endsort"
	if [ `check_done $lastcheck $log` -ne 1 ]; then
		continue
	fi

	# Decompressing; checking whether we're talking about SRA, gzipped FASTQ
	# or normal FASTQ. For the latter two, we make soft links to the target 
	# files so that they appear in the temp (for convenience). If we've got 
	# paired end data, we split the SRA or we find the matching FASTQ and 
	# extract that as well. For the latter, we assume everything is
	# in the same directory, otherwise it'd be harder to code.
	curjob=$curfile":split"
	if [ `check_done $curjob $log` -eq 1 ]; then
		if [[ $curfile =~ \.sra$ ]]; then
			if [ $pet -eq 1 ]; then
	   			fastq-dump -O $temp --split-files $curfile
			else
	   			fastq-dump -O $temp $curfile
	 	fi
		elif [[ $curfile =~ \.(fastq|fq)(\.gz)?$ ]]; then
			if [ $pet -eq 1 ]; then
				matefile=`echo $curfile | sed -r "s/1(\.(fastq|fq)(\.gz)?)$/2\1/"`
				if [[ $matefile == $curfile ]] || [ ! -e $matefile ]; then
					# If the mate doesn't exist, we bail.
 					continue
				fi
				ln $matefile $temp
			fi
			ln $curfile $temp
		else
			echo "Unsupported file format: $curfile"
			exit 1
		fi 

		# FastQC doesn't like the FQ extension, so we will change the name to avoid it.
		for fn in `ls $temp | grep "\.fq"`
		do 
			newloc=`echo $fn | sed "s/\.fq/.fastq/"`
			mv -i $temp/$fn $temp/$newloc
		done
		echo $curjob >> $log
	fi

	# Checking the files. This needs to be done if the file needs *any* processing.
	if [ $pet -eq 1 ]; then
		fastq1=`find $temp | egrep "1\.fastq(\.gz)?$"`
		fastq2=`find $temp | egrep "2\.fastq(\.gz)?$"`
		if [[ $fastq1 == "" ]] || [[ $fastq2 == "" ]]; then
			echo "Paired FastQ files not found"
			exit 1
		fi
		prefix1=`echo $fastq1 | sed -r "s/_[^_]*1\.fastq(\.gz)?$//"`
		prefix2=`echo $fastq2 | sed -r "s/_[^_]*2\.fastq(\.gz)?$//"`
		if [[ $prefix1 != $prefix2 ]]; then
			echo "Prefixes do not match up for the two FastQ files"
			exit 1
		fi
		prefix=`basename $prefix1`
	else
		fastq1=`find $temp | egrep "\.fastq(\.gz)?$"`
		if [[ $fastq1 == "" ]]; then
			echo "FastQ file not found"
			exit 1
		fi
		prefix=`basename $fastq1 | sed -r "s/\.fastq(\.gz)?$//"`
	fi
	iszipped=""
	if [[ $fastq1 =~ \.gz$ ]]; then
		iszipped="--gzFASTQinput"
	fi
	
	# Making QC plots (removing the *.zip file that goes with it).
	curjob=$curfile":plot"
	if [ `check_done $curjob $log` -eq 1 ]; then
		fastqc $fastq1 -o $pics
		if [ $pet -eq 1 ]; then
			fastqc $fastq2 -o $pics
		fi	
		rm $pics/*.zip			  
		echo $curjob >> $log
	fi

	# Aligning files with subread-align. We provide some protection against
	# very short reads, so that the consensus threshold isn't too high.  We peek at
	# the file and check if it's too small.
	curjob=$curfile":align"
	curbam=$temp/$prefix".bam"
	if [ `check_done $curjob $log` -eq 1 ]; then
		conthresh=3
		firstseq=`sed '2q;d' $fastq1`
		if [ ${#firstseq} -lt 45 ]; then
			conthresh=2
		fi

		# Checking for Phred'ness, given that it's not guaranteed to be constant across libraries.
        curphred=`guess_phred $pics $fastq1`
		if [ $curphred -eq 64 ]; then
			curphred=6
		elif [ $curphred -eq 33 ]; then
			curphred=3
		else
			curphred=$phred
		fi
		
		subcmd=/export/share/elvis/bioinf/liao/Subread-1.4.4-p1/bin/subread-align
		if [ $pet -eq 1 ]; then
			$subcmd -i $genome -r $fastq1 -R $fastq2 -o $curbam -P $curphred -m $conthresh --BAMoutput $iszipped
		else
			$subcmd -i $genome -r $fastq1 -o $curbam -P $phred -m $conthresh --BAMoutput $iszipped
		fi
		echo $curjob >> $log
	fi

	# Sorting by mapping coordinate.
	curjob=$curfile":possort"
	altpre=$temp/temp
	altbam=$altpre".bam"
	if [ `check_done $curjob $log` -eq 1 ]; then
		samtools sort $curbam $altpre
		mv $altbam $curbam
		echo $curjob >> $log
	fi

	# Marking duplicate reads.
	curjob=$curfile":dedup"
	if [ `check_done $curjob $log` -eq 1 ]; then
		MarkDuplicates I=$curbam O=$altbam M=$temp/blah.txt TMP_DIR=$vtmp AS=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
		mv $altbam $curbam
		echo $curjob >> $log
	fi

	# Resorting PET data to recover the original order, if requested.
	# Then, moving it to its new home. No point checking for jobness, as
	# we'd have finished up top if that were the case.
	if [ $pet -eq 1 ] && [ $sortbyname -eq 1 ]; then
		samtools sort -n $curbam $bam/$prefix
	else 
		mv $curbam $bam
		samtools index $bam/$prefix".bam"
	fi
	echo $lastcheck >> $log

	# Cleaning out the temporary directories.	
	rm $temp/*
done

###############################################################################
# Erasing the leftover directories, once it's clear that everything has been
# processed successfully (no need to hold onto them).

rmdir $temp
rm -r $vtmp
rm $log   

###############################################################################
###############################################################################
###############################################################################
# Finish.

