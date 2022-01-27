Here is described how to compute a code profile for ECOGEN.

1. Open a terminal in ECOGEN directory.

2. Execute a 'make clean' command.

3. Compile the code with the '-pg' option (see code-profile CXXFLAGS in ECOGEN Makefile).

4. Execute the code either in sequential or in parallel for the test case you want the profile. This will generate a gmon.out file in ECOGEN directory.

5. Install gprof if not already done.

6. For an analysis through a file, execute the following command.
   gprof ECOGEN gmon.out > coverage_and_profile/analysis.txt

7. For a dot-graph visualization, install gprof2dot (https://github.com/jrfonseca/gprof2dot) if not already done and execute a command such as the following.
   gprof ECOGEN gmon.out | gprof2dot | dot -Tpng -o coverage_and_profile/profile_gprof2dot.png

8. Open coverage_and_profile/profile_gprof2dot.png to visualize the profile.