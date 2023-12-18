# Install g++ compiler
conda install -c conda-forge gcc_linux-64 --force-reinstall

# Set environment variables to use the Conda environment's compilers
export CC=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc
export CXX=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++

# compile CNEFinder
cd $SRC_DIR/CNEFinder/
./pre-install.sh
make -f Makefile
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