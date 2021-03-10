#!/bin/bash
#set -e # Exit with nonzero exit code if anything fails

error=0

#Settings
#--------
LIST_TESTS_NONREG="./nonreg/ECOGEN_nonReg_full.list"
#LIST_TESTS_NONREG="./libTests/tests.list"
MAIN_OUTPUT="mainOutput.out"
REPORT_FILE="report.out"
REFERENCE_BRANCH=origin/devel
#REFERENCE_BRANCH="devel"
if [ -z $CI ] 
then
	VALIDATION_BRANCH=$(git symbolic-ref HEAD --short)
else
	VALIDATION_BRANCH=$CI_BUILD_REF_NAME
fi
echo "Branch to validate: $VALIDATION_BRANCH"
validationBranchSimplified=$(echo $VALIDATION_BRANCH | sed -e "s/\//_/g")

#Creating folder for report storage
#----------------------------------
reportFolderName="./nonreg/reports/$(date +%Y%m%d_%H%M)_$validationBranchSimplified"
mkdir -p $reportFolderName

#Cleaning results
#----------------
rm -Rf nonreg/results_reference nonreg/results_validation
echo "nonreg/results_reference and nonreg/results_validation deleted"
mkdir -p results
mv -v results results_save
if [ $? != 0 ]
then
	exit 1
fi

echo "************************************************************************"
echo "--------    RUNNING $VALIDATION_BRANCH BRANCH NONREG TESTS   -----------"
echo "************************************************************************"
echo "Cleaning and compiling..."
make clean > $reportFolderName/output_compile.out
make -j 4 >> $reportFolderName/output_compile.out 2>&1
if [ $? == 0 ]
then
	echo "Compilation OK"
	echo "Running non-regression tests for $VALIDATION_BRANCH branch... More details in main output file ./results/$MAIN_OUTPUT and report file ./results/$REPORT_FILE"
	chmod u+x ./scripts/run.sh
	./scripts/run.sh $LIST_TESTS_NONREG ./results/$MAIN_OUTPUT ./results/$REPORT_FILE
	mv -v results $reportFolderName/results_validation
else
	echo "Compilation failed. Running on $VALIDATION_BRANCH abort."
	error=1
fi

echo "************************************************************************"
echo "----------------   RUNNING $REFERENCE_BRANCH NONREG TESTS   --------------------"
echo "************************************************************************"
echo "Fetching and checking out $REFERENCE_BRANCH branch..."
git fetch > $reportFolderName/output_fetchAndCheckout.out
git branch -D reference >> $reportFolderName/output_fetchAndCheckout.out
git checkout -b reference $REFERENCE_BRANCH >> $reportFolderName/output_fetchAndCheckout.out
if [ $? != 0 ]
then
	echo "Error during checkout on $REFERENCE_BRANCH. Running on $REFERENCE_BRANCH abort."
	error=1
else
	echo "Cleaning and compiling..."
	make clean > $reportFolderName/output_compile_ref.out
	make -j 4 >> $reportFolderName/output_compile_ref.out 2>&1
	if [ $? == 0 ]
	then
		echo "Compilation OK"
		echo "Running non-regression tests for $REFERENCE_BRANCH branch... More details in main output file ./results/$MAIN_OUTPUT and report file ./results/$REPORT_FILE"
		chmod u+x ./scripts/run.sh
		./scripts/run.sh $LIST_TESTS_NONREG ./results/$MAIN_OUTPUT ./results/$REPORT_FILE
		mv -v results $reportFolderName/results_reference
	else
		echo "Compilation failed. Running on $REFERENCE_BRANCH abort."
		error=1
	fi
fi
git checkout $VALIDATION_BRANCH
git branch -D reference >> $reportFolderName/output_fetchAndCheckout.out

echo "************************************************************************"
echo "--------   FINISHING NONREG TESTS (CLEANING/COMPUTING DIFF)  -----------"
echo "************************************************************************"
echo "Cleaning..."
mv -v results_save results

echo "-----------------------------------------"
echo "Computing differences between reference and validation branches..."
diff -x $MAIN_OUTPUT -qr -x 'infoCalcul.out' -x '.DS_Store' $reportFolderName/results_reference $reportFolderName/results_validation > $reportFolderName/diff_nonreg.out || true
if [ -s $reportFolderName/diff_nonreg.out ]
then
	echo "Differences exist: Details available in $reportFolderName/diff_nonreg.out file."
	error=1
else
	echo "Non-regression tests complete: Results present no difference!"
fi

#Sending error code if something went wrong
if [ $error == 1 ]
then
	exit 1
fi
exit 0