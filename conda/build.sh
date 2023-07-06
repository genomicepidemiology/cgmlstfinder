#!/bin/bash

mkdir -p ${PREFIX}/bin

chmod +x cgMLST.py
cp cgMLST.py ${PREFIX}/bin/cgMLST.py

# copy script to download database
chmod +x ${RECIPE_DIR}/download-cgmlstfinder-db.sh
cp ${RECIPE_DIR}/download-cgmlstfinder-db.sh ${PREFIX}/bin/download-cgmlstfinder-db.sh
