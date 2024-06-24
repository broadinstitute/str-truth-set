set -ex

for d in expansion_hunter gangstr hipstr expansion_hunter_denovo  trgt  longtr   straglr  
do
    echo $d
    cd $d;
    make
    cd ..
done
