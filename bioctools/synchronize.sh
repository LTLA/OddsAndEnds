# This will apply changes in repository 1 over to repository 2.

if [[ $# -ne 2 ]]
then
	echo $0 "<source> <destination>"
	exit 
fi

for var in "$@"
do
	if [[ -e $var/.svn ]]
	then
		out=`svn status $var`
	elif [[ -z `cd $var; git rev-parse 2>&1 | grep fatal` ]]
	then
		out=`cd $var; git status $var | egrep "nothing to commit"`
		if [[ -z $out ]] 
		then
			out='failed'
		else 
			out=''
		fi
	else 
		echo "Neither Git or SVN repository detected in $var"
		exit
	fi

	if [[ ! -z $out ]]
	then
		echo "Repository $var is not fully committed"
		exit
	fi
done

rsync -azv --delete --exclude=.svn $1 $2
