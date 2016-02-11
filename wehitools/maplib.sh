###############################################################################
###############################################################################
# This script acts like a BASH library to provide some functions and routines
# for use in mapping. In particular, it loads in the log file and defines 
# a function that checks the log array for jobs that have already been done.
# It also defines a function fo ruse in Phred-checking.

###############################################################################
# This function checks whether the current job (first argument) is present in
# the log file (second argument).

function check_done {
	if [ ! -e $2 ]; then
		echo 1 # Empty log file, so the thing can't be there.
		return 0
	fi
	if grep -qxF $1 $2; then
		echo 0
	else
		echo 1
	fi
	return 0
}

###############################################################################
# This function determines the Phred offset by looking at the FASTQC output, 
# given the picture directory and the name of one of the files.

function guess_phred {
	cleaned=`basename $2 | sed -r "s/\.fastq(\.gz)?$//g"`
	curfile=$1/$cleaned"_fastqc"/fastqc_data.txt
        if [ ! -e $curfile ]; then
		echo "no FASTQC output for $2" >&2
		return 1
        fi
	if grep -q "^Encoding.*Sanger" $curfile; then # Definitely Phred+33
		echo 33
	elif grep -q "^Encoding.*Illumina" $curfile; then # Maybe one or the other, depending on version.
		version=`grep "^Encoding.*Illumina" $curfile | sed "s/.*Illumina[^0-9]*//g"`
		if [ `echo -e "$version\n1.8" | sort -V | head -n1` != "1.8" ]; then
			echo 64
		else
			echo 33
		fi			
	else # Otherwise unknown.
		echo "Can't identify Phred for $2" >&2
		return 1
	fi
	return 0
}

###############################################################################
###############################################################################
# Finish.

