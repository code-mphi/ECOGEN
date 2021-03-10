Here is described how to compute a code coverage for ECOGEN.

1. Open a terminal in ECOGEN/ECOGEN directory.

2. Execute a 'make clean' command.

2. Compile the code with the '-fprofile-arcs -ftest-coverage' options (see code-coverage CXXFLAGS in ECOGEN Makefile).
   This will generate .gcno files in the sources.

3. Execute the code either in sequential or in parallel (recommended) taking care of enabling all the test cases you want to compute the coverage. Note that you can execute the code several times and it will stack up the coverage analysis. For example, enable all the sequential and parallel nonreg test cases and run:
   ./ECOGEN nonreg/ECOGEN_nonreg_sequential.xml && mpirun -n 7 --oversubscribe ./ECOGEN nonreg/ECOGEN_nonreg_parallel.xml
   This will generate .gcda files in the sources.

4. Install lcov if not already done and execute the following command to interprete the gcda files. Note that the "--capture" option may be "--coverage" on some systems.
   lcov --capture --directory . --output-file coverage/coverage.info

5. Execute the following command to generate html files from the previous coverage files.
   genhtml coverage/coverage.info --output-directory coverage/html

6. Open coverage/html/index.html with your browser to have a look at the report.

Note that the .gcno and .gcda files are voluntarily not added in the .gitignore to remind you to delete them after the analysis (and obviously if you do not want to keep them).