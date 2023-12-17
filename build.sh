# copy git directory structure to conda package location
mkdir -p $PREFIX/main $PREFIX/aux $PREFIX/bin $PREFIX/db $PREFIX/lib
cp -r aux/* $PREFIX/aux/
cp -r main/* $PREFIX/main/
cp -r bin/* $PREFIX/bin/
cp -r db/* $PREFIX/db/
cp -r lib/* $PREFIX/lib/

# compile CNEFinder
cd CNEFinder/
./pre-install.sh
make -f Makefile