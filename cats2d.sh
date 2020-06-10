#!/bin/bash


#Copyright (C) 2018 Andreas Mayr
#Licensed under GNU General Public License v3.0 (see http://www.bioinf.jku.at/research/lsc/LICENSE and https://github.com/ml-jku/lsc/blob/master/LICENSE)

rawDataDir=chembl20
sdfFile=smiles.sdf

mkdir mydata/raw/$rawDataDir
mkdir mydata/raw/$rawDataDir/CATS2D
cp smiles.sdf mydata/raw/$rawDataDir
./chemblScript2.sh CATS2D ELEMENT_SYMBOL mydata/raw/$rawDataDir/$sdfFile mydata/raw/$rawDataDir/CATS2D chembl_id STRING_PATTERNS


dirName=mydata/raw/$rawDataDir/CATS2D
outFile=CATS2D_ES.fpf

rm $dirName/../$outFile
for i in `ls $dirName/myout*.csv`; do
  echo $i
  head -n`wc -l $i |cut -f 1 -d " "` $i >> $dirName/../$outFile

done


cp mydata/raw/$rawDataDir/CATS2D/myout1.csv .

