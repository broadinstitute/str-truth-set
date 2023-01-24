set -ex

cd expansion_hunter
make
cd ..
cd gangstr
make
cd ..
cd hipstr
make
cd ..
cd expansion_hunter_denovo
make
cd ..