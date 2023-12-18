# compile CNEFinder
cd $SRC_DIR/CNEFinder/
./pre-install.sh
make -f Makefile CC=$PREFIX/bin/x86_64-conda-linux-gnu-gcc CXX=$PREFIX/bin/x86_64-conda-linux-gnu-g++
mv cnef ../bin/

# copy git directory structure to conda package location
cd $SRC_DIR
mkdir -p $PREFIX/main $PREFIX/aux $PREFIX/bin $PREFIX/db $PREFIX/lib $PREFIX/CNEFinder
cp -r aux/* $PREFIX/aux/
cp -r main/* $PREFIX/main/
cp -r bin/* $PREFIX/bin/
cp -r db/* $PREFIX/db/
cp -r lib/* $PREFIX/lib/
cp -r CNEFinder/* $PREFIX/CNEFinder/