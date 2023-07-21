for ebi_ftp_url_omg in $(cat ../metadata/ebi_ftp_R*.txt ) ;
do
    echo bsub -o $(basename ${ebi_ftp_url_omg}.out) -e $(basename ${ebi_ftp_url_omg}.err) -W 24:00 -J ${ebi_ftp_url_omg} "wget $ebi_ftp_url_omg"
done
