set -u
set -e

# This Bash script is designed to provide support for the Hi-C mapping pipeline
# i.e. to go from *.sra (or FastQ) files to *.bam files.  It does so in a
# manner that is robust to crashes i.e. can be reloaded at each command if
# necessary. I've refrained from putting in an e.g. Python-mediated pipeline as
# I'd have to learn all the coding conventions, and I don't want to do that
# until I know which one to use (or, have one recommended to me).

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

# We need to pull out the location of the script used for mapping. We do
# so and we copy it to the current location, so that re-installs don't 
# affect us.

counterfoil=my_hicMapper.py
alchemist=$temp/$counterfoil
if [ ! -e $alchemist ]; then
	echo "require(diffHic); file.copy(system.file('python', 'presplit_map.py', package='diffHic'), '$alchemist');" | R --no-save
	if [ ! -e $alchemist ]; then
		echo "Extraction of Hi-C mapping script failed."
		exit 1
	fi
fi

# We load some common functions that we need from a BASH mapping library.
# This includes all the logging functions and something to guess the Phred
# offset so that Bowtie2 doesn't throw up.

log=trim.log
source `dirname "${BASH_SOURCE[0]}"`/maplib.sh

###############################################################################
# Now, proceeding to iterate. We make sure we pull out both of the resulting
# *.fastq files from the single *.sra file by using some trickery (i.e. adding
# 1.fastq and 2.fastq to the end of each). Alternatively, we can match up
# Gzipped FASTQ files in each pair. We check whether each one is Phred33 or
# Phred64 using the FASTQC output.
#
# As a side note, some non-duplicate pairs will have the same pair of positions
# if both reads are chimeric and their 5' segments are on the reverse strand, 
# i.e., the recorded positions are the restriction sites. These will not be
# reported as duplicates by MarkDuplicates, as the positions actually used in
# the function are those for the 5' ends of each read. Check out
# < http://sourceforge.net/p/samtools/mailman/samtools-help/thread/4C5C47A0.9060106@eva.mpg.de/ >
# for more details.

for curfile in ${files[@]}; do
	lastcheck=$curfile":endsort"
	if [ `check_done $lastcheck $log` -ne 1 ]; then
		continue
	fi

	# Decompressing.
	curjob=$curfile":split"
	if [ `check_done $curjob $log` -eq 1 ]; then
		compressfind='\.(fastq|fq)\.gz$'
		if [[ $curfile =~ \.sra$ ]]; then
	   		fastq-dump -O $temp --split-files $curfile
		elif [[ $curfile =~ ${compressfind} ]]; then
			if [[ $curfile =~ 1${compressfind} ]]; then
				matefile=`echo $curfile | sed -r "s/1(${compressfind})/2\1/"`
				useme=1
			elif [[ $curfile =~ 2${compressfind} ]]; then
				matefile=`echo $curfile | sed -r "s/2(${compressfind})/1\1/"`
				useme=0
			fi
			# We check that they match up, but we only process when curfile=1.
			if [ ! -e $matefile ]; then
				echo "Mate file for paired-end data not found"
				exit 1
			fi
			if [[ $useme -eq 0 ]]; then
				continue
			fi
			zcat $curfile > $temp/`basename $curfile | sed "s/\.gz$//"`
			zcat $matefile > $temp/`basename $matefile | sed "s/\.gz$//"`
        fi

		# FastQC doesn't like the FQ extension, so we will change the name to avoid it.
		for fn in `ls $temp | grep "\.fq$"`
		do 
			newloc=`echo $fn | sed "s/\.fq$/.fastq/"`
			mv -i $temp/$fn $temp/$newloc
		done
		echo $curjob >> $log
	fi

	# Checking the files. This needs to be done if the file needs *any* processing.
	fastq1=`find $temp | egrep "1\.fastq$"`
	fastq2=`find $temp | egrep "2\.fastq$"`
	if [[ $fastq1 == "" ]] || [[ $fastq2 == "" ]]; then
		echo "Paired FastQ files not found"
		exit 1
	fi
	prefix1=`echo $fastq1 | sed -r "s/_R?1\.fastq$//"`
	prefix2=`echo $fastq2 | sed -r "s/_R?2\.fastq$//"`
	if [[ $prefix1 != $prefix2 ]]; then
		echo "Prefixes do not match up for the two FastQ files"
		exit 1
	fi
	prefix=`basename $prefix1`

	# Making QC plots (removing the *.zip file that goes with it).
	curjob=$curfile":plot"
	if [ `check_done $curjob $log` -eq 1 ]; then
		echo $fastq1 $fastq2
		fastqc $fastq1 -o $pics --extract
		fastqc $fastq2 -o $pics --extract
		rm $pics/*.zip $pics/*.html
		echo $curjob >> $log
	fi

	# Aligning files with the map-and-split approach.
	rawbam=$temp/$prefix".bam"
	curjob=$curfile":align"
	if [ `check_done $curjob $log` -eq 1 ]; then
		# Checking for Phred'ness, given that it's not guaranteed to be constant across libraries.
        	curphred=`guess_phred $pics $fastq1`
		if [ $curphred -eq 0 ]; then
			curphred=$phred
		fi
		python $alchemist -o $rawbam -G $genome -1 $fastq1 -2 $fastq2 --sig $ligsig -P $curphred --cmd "bowtie2 -p 8"
		echo $curjob >> $log
	fi

	# Checking that the raw BAm file is there, named by its prefix.
	if [[ ! -e $rawbam ]]; then
		echo "Raw BAM file not found, name may differ from expected value" 
		exit 1
	fi

	# Fixing mate information.
	curjob=$curfile":fix"
	fixbam=$temp/fixed_$prefix".bam"
	if [ `check_done $curjob $log` -eq 1 ]; then
		FixMateInformation I=$rawbam O=$fixbam TMP_DIR=$vtmp VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
		echo $curjob >> $log
	fi

	# Removing duplicate reads. 
	curjob=$curfile":dedup"
	if [ `check_done $curjob $log` -eq 1 ]; then
		MarkDuplicates I=$fixbam O=$rawbam M=$temp/blah.txt TMP_DIR=$vtmp AS=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
		echo $curjob >> $log
	fi

	# Resorting to recover the original order, and moving it to its new home.
	# No need to check, as we've already checked at the start of this iteration.
	samtools sort -n $rawbam $bam/$prefix
	echo $lastcheck >> $log

	# Cleaning out the temporary directories.	
	for xf in `ls $temp | grep -vxF $counterfoil`; do
		rm $temp/$xf
	done
	rm $vtmp/*
done

###############################################################################
# Erasing the leftover directories, once it's clear that everything has been
# processed successfully (no need to hold onto them).

rm $alchemist
rmdir $temp
rmdir $vtmp
rm $log   

# Adding a ticket to indicate successful completion of the run; and the version numbers involved.

ticket=success.log
if [[ -e $ticket ]]; then
	rm $ticket
fi

set +e
fastqc -v >> $ticket
echo "write(paste('diffHic', packageVersion('diffHic'), 'in', version\$version.string), '$ticket', append=TRUE)" | R --no-save
stored=`cutadapt --version`
printf "Cutadapt version $stored\n" >> $ticket
stored=`bowtie2 --version | head -1 | sed "s/.*version //g"`
printf "Bowtie2 version $stored\n" >> $ticket
stored=`samtools 2>&1 | grep -i "Version:"`
printf "Samtools $stored\n" >> $ticket
stored=`FixMateInformation --version 2>&1`
printf "FixMateInformation version $stored\n" >> $ticket
stored=`MarkDuplicates --version 2>&1`
printf "MarkDuplicates version $stored\n" >> $ticket

###############################################################################
###############################################################################
###############################################################################
# Finish.

