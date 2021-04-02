#!/bin/csh

set i = 0
while ($i < 64)
 qsub build_subgrid_$i.qss
 echo $i
 @ i ++
end

