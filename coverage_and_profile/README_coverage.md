Here is described how to compute a code coverage for ECOGEN.

1. Open a terminal in ECOGEN directory.

2. Execute a 'make clean' command.

3. Compile the code with the '-fprofile-arcs -ftest-coverage' options (see code-coverage CXXFLAGS in ECOGEN Makefile).
   This will generate .gcno files in the sources.

4. Execute the code either in sequential or in parallel (recommended) taking care of enabling all the test cases you want to compute the coverage. Note that you can execute the code several times and it will stack up the coverage analysis. For example, run the following for all the sequential and parallel nonreg test cases:
   ./scripts/run.sh ./nonreg/ECOGEN_nonReg_full.list
   This will generate .gcda files in the sources.

5. Install lcov if not already done and execute the following command to interprete the gcda files. Note that the "--capture" option may be "--coverage" on some systems. You may also add the following options "--ignore-errors inconsistent,inconsistent,gcov,gcov --filter range".
   lcov --capture --directory . --output-file coverage_and_profile/coverage.info

6. Execute the following command to generate html files from the previous coverage files. You may also add the following options "--ignore-errors inconsistent".
   genhtml coverage_and_profile/coverage.info --output-directory coverage_and_profile/html

7. Open coverage_and_profile/html/index.html with your browser to have a look at the report.

Note that the .gcno and .gcda files are voluntarily not added in the .gitignore to remind you to delete them after the analysis (and obviously if you do not want to keep them).