#!/bin/sh

test_name='test_001'

testpath=`pwd`
fortran_src="https://raw.githubusercontent.com/gschnabel/DATP-Fortran/1bc1bd87d7fef967c8718fe88417320e98cae908/DATP.FOR"
python_src="$testpath/../../source/DATP.py"
numdiff_exe="${NUMDIFF_EXE-numdiff}"

# check if numdiff exiss
which $numdiff_exe > /dev/null
if [ $? -ne 0 ]; then
    echo 'Cannot find numdiff program which is required for the numerical comparison.'
    echo 'You can download it from http://www.nongnu.org/numdiff/'
    exit 1
fi

# change into test directory
cd $testpath
retval=$?
if [ $retval -ne 0 ]; then
    echo Cannot change into test directory
    exit 1
fi

# check if test in correctly named directory
echo $testpath | grep "${test_name}$" > /dev/null
retval=$?
if [ $retval -ne 0 ]; then
    echo This test is run in the wrong directory
    exit 1
fi

# clean previous results
rm -rf result

# create result
mkdir result
cd result

mkdir fortran
cd fortran
echo retrieving Fortran DATP code...
wget --quiet $fortran_src
if [ $? -ne 0 ]; then
    echo Could not download Fortran DATP code
    exit 1
fi

echo generating Fortran result...
# change a few format descriptors
sed -i \
    -e 's!^\( *100 *FORMAT(\) *[^)]*\().*$\)!\1 2E14.6\2!' \
    -e 's!^\( *200 *FORMAT(\) *[^)]*\().*$\)!\1 2E14.6,12F12.7\2!' \
    -e 's!^\( *290 *FORMAT(\) *[^)]*\().*$\)!\1 2E14.6,12F12.7,F9.5\2!' \
    DATP.FOR
# and compile
gfortran -o DATP-gfort DATP.FOR 
cp $testpath/input/* .
./DATP-gfort > output

echo generating Python result...
cd $testpath/result
mkdir python
cd python
cp $testpath/input/* .
export TEST_DATP=yes
python $python_src > output

# compare the output files
cd $testpath/result
$numdiff_exe -r 5e-5 fortran/DAT.RES python/DAT.RES > numdiff.out

retval=$?
if [ $retval -eq 0 ]; then
    echo test successful
else
    echo test failed
fi

exit $retval
