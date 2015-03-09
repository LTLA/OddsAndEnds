# This script checks whether or not the version number was updated. 

if [[ $# -ne 1 ]]
then
	echo `basename $0` "<repository>"
	exit 
fi

if [[ -e $1/.svn ]]
then

	#echo -n "SVN password: "
	#read -s password
	#echo
		
	file=`echo $1 | sed "s/\/$//"`
	if [[ -h $file ]]
	then
		file=`readlink $file`
	fi

	# Updating.
	#extra="--password $password"
	extra=""
	svn up $file $extra > /dev/null

	# Comparing changes on Version: in description, 
	change1=`svn info $file | grep "Last Changed Rev" | sed "s/.*: //"`
	change2=`svn annotate $file/DESCRIPTION $extra | grep "Version:" | sed -E "s/^\s+//g" | cut -f1 -d " "`
elif [[ -z `cd $1; git rev-parse 2>&1 | grep fatal` ]]
then
	cd $1
	change1=`git log -n 1 | grep "commit" | cut -f2 -d ' '`
	change2=`git blame -l DESCRIPTION | grep "Version" | cut -f1 -d ' '`
else 
	echo "Neither Git or SVN repository detected"
	exit
fi

if [[ $change1 != $change2 ]]
then
	echo "Old version number ($change2) in directory '$1' ($change1)"
fi


