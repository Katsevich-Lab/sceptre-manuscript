#!/bin/bash

# This script builds a package in the code directory. It takes as arguments (1) the name of the package, and (2) the name of the machine.

# save the arguments in variables
package=$1
machine=$2

# obtain the code directory
code_dir=$(bash get_file_paths.bash $machine code)

# change directories to the code directory
if [ $package = sceptre ]
then
  to_cd=$code_dir"/sceptre_package"
fi
if [ $package = katsevich2020 ]
then
  to_cd=$code_dir"/sceptre_paper"
fi

cd $to_cd

# Obtain the version of the package through an awk call on DESCRIPTION
version=$(awk '/Version/ {print $2}' $package"/"DESCRIPTION)

# build the package
R CMD build $package

# install the package
R CMD INSTALL $package"_"$version.tar.gz

# remove the tar.gz file
rm $package"_"$version".tar.gz"

# go back to utilities
cd $code_dir"/sceptre_paper/utilities"
