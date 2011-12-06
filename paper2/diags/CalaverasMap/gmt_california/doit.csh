#!/bin/csh

gmtset PAPER_MEDIA A4

set fig1 = "-Xa1i -Ya1i"
set bounds1 = "-R-124/-120/36/40"
set proj1 = "-JM5i"
set anot1 = " -B1g30"


#set bounds2 = "-R-121.8/-121.5/37.1/37.4"
set bounds2 = "-R-600/600/-600/800"
set proj2 = "-JX1.2i/1.2i"
set anot2 = " -BnSeW400"
set fig2 = "-Xa1.6i -Ya1.5i"


# plot topo
grdraster 1 $bounds1 $proj1  -Gttopo2m.grd

# resample the seafloor
grd2xyz ttopo2m.grd>! ttopo2m.xyz
surface ttopo2m.xyz $bounds1 -Gttopo.grd -I1m -T1
# calculate the gradient
grdgradient ttopo.grd -V -A90 -Gtemp.grd
# illuminate the grid
grdmath temp.grd 9000. / = gtopo_grad.grd
\rm temp.grd
# plot the seafloor
grdimage $fig1  ttopo.grd -Igtopo_grad.grd $anot1 $bounds1 $proj1 -Cgrltoposea.cpt  -K -P -V >! tmp.ps

# clip the dry areas

grdraster 22 -I0.5m $bounds1  -Gttopo30.grd

pscoast $fig1 $bounds1 $proj1 -O -K -Dh -A20 -Gc -P >> tmp.ps

# start to plot the dem for the dry areas
# calcualte the gradient
grdgradient ttopo30.grd -V -A90  -Gtemp.grd
# illuminate the grid
grdmath temp.grd 9000. / = gtopo_grad.grd
# plot the land
grdimage $fig1 ttopo30.grd -Igtopo_grad.grd $bounds1 $proj1 -Cgrltopoland.cpt -O -K -P >> tmp.ps
\rm temp.grd

# plot the coastline
pscoast $fig1 $bounds1 $proj1  -W3 -A20 -Dh -O -P -K >> tmp.ps

#switch off the clipping of the dry areas
pscoast $fig1 -Q -O -K $bounds1 $proj1>> tmp.ps



#plotting the stations

# plot stations
awk '{ print $3 "  " $2}' station.dat   | psxy $fig1  $proj1 $bounds1   -St0.05i -O  -P -K -Gred >> tmp.ps



# create the white filled box within which to place in inset
echo "236 36" > poly.txt
echo "237.55 36" >> poly.txt
echo "237.55 37.2" >> poly.txt
echo "236 37.2" >> poly.txt
psxy $fig1 $bounds1 $proj1 poly.txt -L -O -P -K -Gwhite>> tmp.ps


# join insert to area specified on large map
echo "237.55     37.2" > poly2.txt
echo "-121.6620  37.2847" >> poly2.txt
echo "-121.6620  37.2945" >> poly2.txt
echo "237.55     36" >> poly2.txt
psxy $fig1 $bounds1 $proj1 poly2.txt -O -P -K  >> tmp.ps





# COLOR_MODEL = RGB
echo "0 215 025 028  1 215 025 028" > mag.cpt
echo "1 253 174 097  2 253 174 097" >> mag.cpt
echo "2 255 255 191  3 255 255 191" >> mag.cpt
echo "3 171 221 164  4 171 221 164" >> mag.cpt
echo "4 043 131 186  10 043 131 186 " >> mag.cpt
# plot the events by magnitude
awk '{ print $5 "  " $6 "  " $17}' hypoDD.reloc | psxy $fig2 $anot2  $bounds2 $proj2 -Sc0.01i -O -P -Cmag.cpt >> tmp.ps

# Plot a polygon surrounding the events to be used. 
echo "-577.4390    550.6190 " > poly3.txt
echo " 31.0976    -477.7166 " >> poly3.txt
echo "102.1098    -406.7044 " >> poly3.txt
echo "-506.4268    621.6312 " >>poly3.txt
echo "-577.4390    550.6190 " >> poly3.txt
awk '{ print $1 "  " $2}' poly3.txt | psxy $fig2 $bounds2 $proj2 -P >> tmp.ps   #-L -O -P -K  >> tmp.ps

# plot the colourmap
#psscale -Cdepth.cpt -Ba10f5::/:M: -D15/3/6/0.6 -O  -P >> tmp.ps

ps2epsi tmp.ps
gv tmp.epsi
