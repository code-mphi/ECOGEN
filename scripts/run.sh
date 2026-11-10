#!/bin/bash

# -------------------------------------------------------------------------
#
# Print help command line options.
#
# -------------------------------------------------------------------------
usage() {
  echo ""
  echo "Usage: $(basename "$0") [-e ecogen_exec] [-m mail_address] [-l test_list] [-o output_file] [-r report_file]"
  echo ""
  echo "    Run tests listed in <test_list> file, save exit statuses of each test in <report_file>"\
       "file and test terminal outputs in <output_file> file."
  echo "    If a <mail_address> is provided, send <report_file> by email."
  echo ""
  echo "        -e <ecogen_exec> path and name of ECOGEN executable."\
       "Defaults to './ECOGEN'."
  echo "        -l <test_lists> path and name of file with test list."\
       "Defaults to './libTests/tests.list'."
  echo "        -m <mail_address> recipient to send final report to. Email is not sent if not provided."
  echo "        -o <ouput_file> path and name of file in which are redirected test terminal outputs."\
       "Defaults to './results/mainOutput.out'."
  echo "        -r <report_file> path and name of file in which are save test exit statuses."\
       "Defaults to './results/report.out'."
  echo ""
}



# Parse command line options
while getopts ':he:l:m:o:r:' OPTION; do
  case "$OPTION" in
    h)
      usage
      exit 0
      ;;
    e)
      ecogenExec=$OPTARG
      ;;
    l)
      listTestFileName=$OPTARG
      ;;
    m)
      mail=$OPTARG
      ;;
    o)
      mainOutput=$OPTARG
      ;;
    r)
      reportFile=$OPTARG
      ;;
    \?)
      echo ""
      echo "Error: Illegal option -$OPTARG"
      usage
      exit 1
      ;;
    :)
      echo ""
      echo "Error: Option -$OPTARG requires an argument"
      usage
      ;;
  esac
done

# Assign default values if parameters are unset
listTestFileName=${listTestFileName:="./libTests/tests.list"}
mainOutput=${mainOutput:="./results/mainOutput.out"}
reportFile=${reportFile:="./results/report.out"}
mail=${mail:="nomail"}
ecogenExec=${ecogenExec:="./ECOGEN"}

# Test existence of ecogen executable
if [ ! -f "${ecogenExec}" ]; then
  echo "Error: ${ecogenExec} executable not found."
  exit 1
fi

# Test existence of test list
if [ ! -f "${listTestFileName}" ]; then
  echo "Error: ${listTestFileName} file not found."
  exit 1
fi

#Saving main input file ECOGEN.xml
#---------------------------------
fich="ECOGEN.xml"
mv ./"$fich" ./mainInput.save
mkdir -p results

#Reading and storing test-case list
#----------------------------------
nb_tests=$(awk '!/^#/ && NF > 0 { count++ } END { print count+0 }' "$listTestFileName")

while IFS=' ' read -r c1 c2 c3 || [[ -n "$c1" ]]; do
  # Skip comments (line starting with '#') and blank lines
  if [[ "$c1" =~ ^# ]] || [[ -z "$c1" ]]; then
    continue
  fi

  if [[ -n "$c1" && -n "$c2" ]]
  then
    testsList+=("$c1")  # Catching test folders
    ncpu+=("$c2")       # Catching core numbers
    retval+=("${c3:-0}") # Expected return value (filled by 0 if not provided)
  fi
done < "$listTestFileName"

rm -f "$mainOutput"
touch "$mainOutput"
rm -f "$reportFile"
touch "$reportFile"
# file to save the list of failed tests, so that only the latter can be rerun
failedList="./results/rerunFailed.list"
rm -f "$failedList"
touch "$failedList"

#Set Gmsh options according to its version
#-----------------------------------------
#By default Gmsh version >4 produces mesh file
#not compatible with ECOGEN. Specific options are set
# to make the mesh file compliant.
gmshOptions=""
gmshOptionsV2=""
gmshOptionsV4=""
if gmsh --version &> "gmsh.txt"
then
  gmshVersion=$(cut -d"." -f1 < gmsh.txt)
  if [ "$gmshVersion" -ge "4" ]
  then
    gmshOptionsV2="-part_ghosts -format msh22"
    gmshOptionsV4="-part_ghosts"
  fi
fi
rm "gmsh.txt"

#Looping on test cases
#---------------------
exit_status=0
testNum=0
for ((k=0;k<"$nb_tests";k++)) #array begin to 0
do
  testNum=$((testNum + 1))
  status="succeed"

  printf 'Test case %s running on %s cores...' "${testsList[$k]}" "${ncpu[$k]}"

  #Running gmsh if necessary
  #*************************
  meshType=$(sed -n 's/<type structure="//p' "${testsList[$k]}mesh.xml" | sed -n 's/"\/>//p' | tr -d '[:space:]')
  meshType=$(echo "$meshType" | tr '[:upper:]' '[:lower:]') # Convert to lowercase
  if [ "$meshType" == "unstructured" ]
  then
    mshFile=$(sed 's/\t//g' "${testsList[$k]}mesh.xml" | sed 's/ //g' | sed -n 's/^<filename="//p' | sed 's/\".*//')
    geoFile=$(echo "$mshFile" | sed 's/msh/geo/')

    # For 'mshfileEncoding/.*V4' test cases we want to generate meshes at V4 gmsh file format if available
    case "${testsList[$k]}" in
      *"mshfileEncoding/"*"V4"*)
        # Options for V4 format if available
        gmshOptions="$gmshOptionsV4"
        ;;
      *)
        # Options for V2 format
        gmshOptions="$gmshOptionsV2"
        ;;
    esac

    if ! eval gmsh -3 -part "${ncpu[$k]}" "$geoFile" "$gmshOptions" >> "$mainOutput" 2>&1
    then
      status="failed: mesh generation failed"
    fi

    # Special treatment to test file encoding
    case "${testsList[$k]}" in
      *"mshfileEncoding/dos"*)
        # force 'dos' encoding by inserting CR (carriage return) before LF (line feed)
        sed -i'' 's/$/\r/' "$mshFile" >> "$mainOutput" 2>&1
        ;;
      *"mshfileEncoding/oldMac"*)
        # force old 'mac' encoding by replacing LF (line feed) by CR (carriage return)
        tr '\n' '\r' < "$mshFile" > "$mshFile.tmp" && mv  "$mshFile.tmp" "$mshFile"  >> "$mainOutput" 2>&1
        ;;
    esac
  fi

  if [ "$status" == "succeed" ]
  then
    #Writing main input file ECOGEN.xml
    #**********************************
    echo "<?xml version='1.0' encoding='UTF-8' standalone='yes'?>" > "$fich"
    echo "<ecogen>" >> "$fich"
    chaine="  <testCase>${testsList[$k]}</testCase>"
    echo "$chaine" >> "$fich"
    echo "</ecogen>" >> "$fich"

    #Running ECOGEN
    #**************
    mpirun -np "${ncpu[$k]}" "${ecogenExec}" >> "$mainOutput" 2>&1

    #Writting report
    #***************
    errorCode="$?"
    if [ "$errorCode" != "${retval[$k]}" ]
    then
      status="failed: error $errorCode - expected: ${retval[$k]}"
      exit_status=1
      echo "${testsList[$k]} ${ncpu[$k]} ${retval[$k]}" >> "$failedList"
    else
      status="succeed"
    fi
  fi

  echo "Test case: '${testsList[$k]}' $status" >> "$reportFile"
  printf '  %s\n' "$status"

done

#Reloading main input file ECOGEN.xml
#------------------------------------
rm -f "$fich"
mv ./mainInput.save ./"$fich"

#Ending and sending end e-mail
#-----------------------------
sujet="Job_Complete"
message="Hi Master. My job is complete."
echo "$message"
if [ "$mail" != "nomail" ]
then
  (echo "$message"; cat "$reportFile") | mail -s "$sujet" "$mail"
  echo "Mail sent to: $mail"
fi

exit "$exit_status"
