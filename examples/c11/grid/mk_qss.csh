#!/bin/csh

set i = 0
while ($i < 64)
 cat template_build_subgrid.qss | sed "s/NNN/$i/" > build_subgrid_$i.qss
 @ i ++
end

