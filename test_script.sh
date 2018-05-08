#!/bin/bash
for number in {1..1217}
do
echo $number
./vari_find_paths $number 90000 500 31.reads.dbg restriction_nodes1en31pe
done
exit 0
