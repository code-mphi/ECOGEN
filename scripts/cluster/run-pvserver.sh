#!/bin/bash

# Start a ParaView server
# For more information see tutorial 'ParaView server'

start()
{	
	if [ -f $path_save ]
	then 
		# Read ParaView path
		path_paraview=$(cat "$path_save")

		# Set path to executables
		cmd_mpiexec=$path_paraview"mpiexec"
		cmd_pvserver=$path_paraview"pvserver"
	
		# Command
		$cmd_mpiexec -np $1 $cmd_pvserver --force-offscreen-rendering

	else
		echo 'ParaView path is not set, could not start server.'
	fi
}

set_path()
{
	# Path is stored into .path file in current directory.
	# Be aware that depending on location of this script sudo 
	# permissions might be required.
	
	# Check/edit if ending of given path is '/'
	if [ "${1: -1}" != "/" ]
	then 
		path_to_set=$1"/"
	else
		path_to_set=$1
	fi

	echo "$path_to_set" > "$path_save"
	if [ $? == 0 ]
	then
		echo "ParaView path has been set to $1."
	else
		echo 'Could not set path (sudo access required).'
	fi
}

help()
{
	# Display help
	echo 'Start a ParaView server. Requires to set path to ParaView bin folder.'
	echo
	echo 'Usage: run-pvserver [OPTION]'
	echo 'Options:'
	echo '  -h,  --help           Print help.'
	echo '  -np, --number-cores   Select the number of cores (default is 1 when not provided).'
	echo '  -p,  --path           Set path to ParaView executables'
	echo
	echo 'Example to set path and run server for 1st time:'
	echo '  $ ./run-pvserver.sh --path /opt/ParaView-5.8.0-MPI-Linux-Python3.7-64bit/bin/'
	echo '  Path has been set to /opt/ParaView-5.8.0-MPI-Linux-Python3.7-64bit/bin/.'
	echo '  $ ./run-pvserver # Using default number of core, i.e. single core'
	echo '  Waiting for client...'
	echo '  Connection URL: cs://cid:11111'
	echo '  Accepting connection(s): cid:11111'
}

# Set default number of cores
core=1

# Current path of the script 
path_save=$(dirname "$0")/.path

# Main
if [ $# -gt 0 ]
then
	case "$1" in
	-np | --number-cores)
		if [ $# -eq 2 ]
		then
			core=$2
		fi
		start $core 
		exit 0
		;;
	-p | --path)
		if [ $# -eq 2 ]
		then
			set_path $2
		else
			echo 'Path could not be set.'
		fi
		exit
		;;
	-h | --help)
		help
		exit 0
		;;
esac
else
	start $core
	exit 0
fi
