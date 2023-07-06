#!/usr/bin/env bash

echo "Downloading lastest version of the cgMLSTFinder database to current directory..."

mkdir cgmlstfinder_db
cd cgmlstfinder_db

wget https://bitbucket.org/genomicepidemiology/cgmlstfinder_db/get/master.tar.gz
tar -xvf master.tar.gz --strip-components 1

echo "Installing the cgMLSTFinder database with KMA"
python INSTALL.py

echo "The cgMLSTFinder database has been downloaded and installed."

exit 0