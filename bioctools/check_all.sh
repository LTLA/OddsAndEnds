##################################################
# Checking options:

useVG=0
RCMD=R

OPTIND=1
while getopts ":vr:" opt; do
	case "$opt" in
		v)
			useVG=1
			;;
		r)
			RCMD=$OPTARG
			;;
	esac
done
shift "$((OPTIND-1))"

if [[ $# -ne 1 ]]
then
	echo `basename $0` "[-v] [-r RCMD] <repository>"
	exit
fi
repo=$1

##################################################
# Making a new directory in which to store results.

set -e
set -u
packname=`cat $repo/DESCRIPTION | grep "Package" | sed "s/Package: //"`
vnum=`cat $repo/DESCRIPTION | grep "Version" | sed "s/Version: //"`
tarball=${packname}_${vnum}.tar.gz
curdir=`pwd`

if [[ ! -e "Checks" ]]
then
	mkdir Checks
fi

checkdir=Checks/$repo
if [[ $useVG -eq 1 ]]
then
	checkdir=${checkdir%/}-valgrind
fi
if [[ !  -e $checkdir ]]
then
	mkdir $checkdir
fi

cd $checkdir

##################################################
# Running the check results.

logfile=`pwd`/check.log
if [[ -e $logfile ]]
then
	rm $logfile
fi
if [[ -n `find -maxdepth 1 | grep "${packname}.*.tar.gz"` ]]
then
	rm ${packname}_*.tar.gz
fi
$RCMD CMD build $curdir/$repo >> $logfile 2>&1
echo >> $logfile

# Running CHECK.
START=`date +%s%N`
if [[ $useVG -eq 1 ]]
then
	$RCMD CMD check --use-valgrind $tarball >> $logfile 2>&1
else
	$RCMD CMD check $tarball >> $logfile 2>&1
fi
END=`date +%s%N`
ELAPSED=`echo "scale=8; ($END - $START) / 1000000000" | bc`
echo "Elapsed time for CHECK: $ELAPSED seconds" >> $logfile
echo >> $logfile

# Installing with warnings attached, if possible.
locale=`echo 'cat(paste0(R.home(), "\n"))' | $RCMD --no-save --vanilla --silent | grep -v ">"`
confile=$locale/etc/Makeconf
newconf=${confile}.temp
if [ -w $confile ]
then
	if [ ! -e $newconf ]
	then
		# If we got interrupted last time before we could do the final move, we don't want to overwrite
		# the '.temp' file with our new settings. So, we don't copy if '.temp' already exists.
		cp $confile $newconf
	fi
	cat $newconf | sed "s/^CXXFLAGS.*/CXXFLAGS = -g -O0 \$(LTO) -Wall -Wextra -pedantic/" > $confile
	$RCMD CMD INSTALL $tarball >> $logfile 2>&1
	mv $newconf $confile
else
	$RCMD CMD INSTALL $tarball >> $logfile 2>&1
fi

# Running custom tests.
cd ${packname}.Rcheck/${packname}/tests
if [[ $useVG -eq 1 ]]
then
	bash run_all_tests.sh $RCMD -d valgrind >> $logfile
else
	bash run_all_tests.sh $RCMD >> $logfile
fi
echo >> $logfile
cd ../../..

# Running BiocCheck.
$RCMD CMD BiocCheck $tarball >> $logfile 2>&1
echo >> $logfile

