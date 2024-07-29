set -x

for d in trgt  longtr   straglr   vamos  expansion_hunter gangstr hipstr expansion_hunter_denovo  
do
    echo $d
    cd $d;
    make
    cd ..
done
