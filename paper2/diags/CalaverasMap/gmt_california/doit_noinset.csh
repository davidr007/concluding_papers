#!/bin/csh

gmtset PAPER_MEDIA A4

set bounds = "-R-124/-120/36/40"
set proj = "-JM5i"
set anot = " -BNESWa1f0.5/a1f0.5"
set anot = " -B1g30"

#set fig1 = "-Xa1i -Ya3i"
#set fig2 = "-Xa1i -Ya1i"


# plot topo
grdraster 1 $bounds $proj  -Gttopo2m.grd

# resample the seafloor
grd2xyz ttopo2m.grd>! ttopo2m.xyz
surface ttopo2m.xyz $bounds -Gttopo.grd -I1m -T1
# calculate the gradient
grdgradient ttopo.grd -V -A90 -Gtemp.grd
# illuminate the grid
grdmath temp.grd 9000. / = gtopo_grad.grd
\rm temp.grd
# plot the seafloor
grdimage  ttopo.grd -Igtopo_grad.grd $anot $bounds $proj -Cgrltoposea.cpt  -K -P -V >! tmp.ps

# clip the dry areas

grdraster 22 -I0.5m $bounds  -Gttopo30.grd

pscoast $bounds $proj -O -K -Dh -A20 -Gc -P >> tmp.ps

# start to plot the dem for the dry areas
# calcualte the gradient
grdgradient ttopo30.grd -V -A90  -Gtemp.grd
# illuminate the grid
grdmath temp.grd 9000. / = gtopo_grad.grd
# plot the land
grdimage ttopo30.grd -Igtopo_grad.grd $bounds $proj -Cgrltopoland.cpt -O -K -P >> tmp.ps
\rm temp.grd

# plot the coastline
pscoast $bounds $proj  -W3 -A20 -Dh -O -P -K >> tmp.ps

#switch off the clipping of the dry areas
pscoast -Q -O -K $bounds $proj>> tmp.ps



#plotting the stations

# plot stations
awk '{ print $3 "  " $2}' station.dat   | psxy  $proj $bounds   -St0.05i -O  -P -K -Gred >> tmp.ps


# COLOR_MODEL = RGB
echo "0 215 025 028  1 215 025 028" > mag.cpt
echo "1 253 174 097  2 253 174 097" >> mag.cpt
echo "2 255 255 191  3 255 255 191" >> mag.cpt
echo "3 171 221 164  4 171 221 164" >> mag.cpt
echo "4 043 131 186  10 043 131 186 " >> mag.cpt

#-K tells GMT that another plotting command follows...
# plot the points
awk '{ print $3 "  " $2 "  " $8}' hypoDD.reloc | psxy  $bounds $proj -Sc0.01i -O -P -Cmag.cpt -O -P>> tmp.ps

# plot the colourmap
#psscale -Cdepth.cpt -Ba10f5::/:M: -D15/3/6/0.6 -O  -P >> tmp.ps


ps2epsi tmp.ps
gv tmp.epsi
