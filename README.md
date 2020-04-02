## DATP - Python

The [DATP](https://github.com/gschnabel/DATP-Fortran) code preprocesses 
experimental data sets prior to fitting with the GMAP code.
Both codes have been developed by Wolfgang P. Poenitz in Fortran 77.
As Fortran 77 is not that much used anymore nowadays, the conversion
to Python will facilitate future developments.
Due to the large ecosystem of extension modules around Python, the code
logic of DATP and GMAP can be formulated on a higher level, thereby
reducing code complexity. Moreover, also interfacing with other
applications, such as web applications and databases will become easier.

### Approach to conversion
Because the code is used for the evaluation of neutron cross section
standards, it is very important that the Python version is functionally
equivalent to the Fortran version. Only then can methodological
improvements in the new code be made and their effect on the results
studied.

For this reason, in a first step the code was translated on a line by
line basis to Python keeping exactly the same logic and level of
abstraction as present in the Fortran code, even keeping `goto`
statements which had to be introduced into the Python code using
a bytecode hack implemented by the [goto-statement] module. 
Also the [fortranformat] module proved invaluably useful due to
heavy input/output operations and associated FORMAT edit descriptors
which do not have a perfect counterpart in Python.
 
To ensure the functional equivalence of the Python code and the
Fortran code, a test case was written that compares the numerical
output of both codes. As Python works natively with double
precision floats whereas the Fortran code with single precision
and, in addition, Python and Fortran implement different rounding
schemes, small numerical differences in the outputs have to
be expected. The test case therefore relies on the [numdiff]
tool to compare the outputs which ignores small insignificant 
differences. In the test it is set up to accept relative difference
below 5x10^5.

Assuming that you cloned this repository, in order to run the
test from a command line in Linux, you have to change into the
directory `tests/test_001`. Being there execute:
```
export NUMDIFF_EXE=<path to numdiff binary>
./run_test.sh
```
The shell script will return a text message on the console
whether outputs are equivalent and also indicate the result
with the return code, i.e., 0 for equivalent output and
1 if differences were detected. In the latter case, the file
in `/result/numdiff.out` contains details on the location of
the differences in `result/fortran/DAT.RES` and
`result/python/DAT.RES`.

[goto-statement]: https://pypi.org/project/goto-statement/
[fortranformat]: https://pypi.org/project/fortranformat/ 
[numdiff]: http://www.nongnu.org/numdiff/

### Future developments
Now that a first draft version of the Python code is available,
with large probability fully functionally equivalent to the
Fortran version due to implementing exactly the same program flow
and the equivalence of the output being tested, the Python
will be incrementally refactored, always validating the effect
of updates against the output of the original Fortran version.
