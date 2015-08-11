#!/bin/bash

name="rescan-052k"
data="/home/pixie16/isoldeData/is599/052k/"
confs="/home/pixie16/svp/OptionalConfigurations/is599"

make && rm -f $name.*

ln -s -f $confs/Config.xml.is599_52k_00 Config.xml
cmd="zero\nban 051k.ban\nfile $data/is599_52k_00.ldf\ngo\nend"
echo -e $cmd | ./pixie_ldf_c $name

ln -s -f $confs/Config.xml.is599.52k_01 Config.xml
cmd="ban 051k.ban\nfile $data/is599_52k_01.ldf\ngo\nend"
echo -e $cmd | ./pixie_ldf_c $name

cmd="ban 051k.ban\nfile $data/is599_52k_02.ldf\ngo\nend"
echo -e $cmd | ./pixie_ldf_c $name

cmd="ban 051k.ban\nfile $data/is599_52k_0.ldf\ngo\nend"
echo -e $cmd | ./pixie_ldf_c $name