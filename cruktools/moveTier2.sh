if [[ $1 != "from" ]] && [[ $1 != "to" ]]
then
    echo "Need from/to specification"
    exit 1
fi

if [[ ${#2} -eq 0 ]]
then
    echo "Need directory to get or put stuff on"
    exit 1
fi

if [[ ${#3} -eq 0 ]]
then
    echo "Need mask to determine what to transfer"
    exit 1
fi

if [[ $1 == "from" ]]
then

# Use here documents to feed in commands.
smbclient //jmlab-data/jmlab -U lun01 << SMBCOMMANDS
prompt
cd $2
mget $3
SMBCOMMANDS

elif [[ $1 == "to" ]]
then

smbclient //jmlab-data/jmlab -U lun01 << SMBCOMMANDS
prompt
cd $2
mput $3
SMBCOMMANDS

fi

exit 0
