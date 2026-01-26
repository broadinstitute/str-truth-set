set -x

for d in hipstr longtr  # expansion_hunter trgt  # gangstr expansion_hunter_denovo  
do
    echo $d
    cd $d;
    make
    cd ..
done
