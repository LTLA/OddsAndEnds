for x in $@
do 
    firstthree=$(echo $x | cut -c 1-6)
    ascp -i /usr/local/bioinfsoftware/aspera/current/etc/asperaweb_id_dsa.openssh -k1 -Tr -l100m anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/${firstthree}/$x .
    mv $x/* .
    rmdir $x
done
