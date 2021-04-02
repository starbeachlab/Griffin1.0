#!/bin/bash

tmp="counter.tmp";

value=0;

if [ -e $tmp ]
then
    value=$(cat $tmp | awk '{print $1+1}')
fi

####  ADJUST HERE THE FREQUENCY OF EXTRACTING FORCES FOR VISUALIZATION WITH VMD  ###
j=$(expr $value % 1000)
#####  END  #####

if [ "$j" == "0" ]; then 

####### ADJUST HERE FILENAMES AND LINES TO EXTRACT, e.g. the molecules for which forces shall be visualized  ##########
perl /home/rene/workspace/perl/extract_lines.pl namdx.tmp namdf.tmp forces.vmd 1475 1608 8443 8576
#####  END  #####

fi

#cp namdf.tmp namdforces.$value.txt
#cp namdx.tmp namdcoor.$value.txt

#### ADJUST PATH TO MESSENGER ####
GRIFFIN/bin/griffin_messenger.0.34.exe -forces

echo $value > $tmp

