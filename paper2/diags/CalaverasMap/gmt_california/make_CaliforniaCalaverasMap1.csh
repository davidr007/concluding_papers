#!/bin/csh

gmtset PAPER_MEDIA A4
gmtset TICK_LENGTH 0.1c

set fig1 = "-Xa1i -Ya1i"
set bounds1 = "-R-124/-120/36/40"
set proj1 = "-JM5i"
set anot1 = " -B1g30"


#set bounds2 = "-R-121.8/-121.5/37.1/37.4"
set bounds2 = "-R-600/600/-600/800"
set proj2 = "-JX1.4i/1.4i"
set anot2 = " -BnSeW400"
set fig2 = "-Xa1.4i -Ya1.3i"


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
#grd2cpt ttopo.grd -Cgray -L0/10000 -S0/200/20 >mydata.cpt
#grd2cpt ttopo.grd -Cgray -L0/10000 -S0/2/50 >mydata.cpt
grd2cpt ttopo.grd -Cgray  -S-10/1000/100000 >mydata.cpt
grdimage $fig1  ttopo.grd -Igtopo_grad.grd $anot1 $bounds1 $proj1 -Cgrltoposea.cpt  -K -P -V >! CaliforniaCalaverasMap1.ps

# clip the dry areas

grdraster 22 -I0.5m $bounds1  -Gttopo30.grd

pscoast $fig1 $bounds1 $proj1 -O -K -Dh -A20 -Gc -P  >> CaliforniaCalaverasMap1.ps
#pscoast $fig1 $bounds1 $proj1 -O -K -Dh -A20 -Gc -P -Lf-125/37.3/37.3/150k >> CaliforniaCalaverasMap1.ps
#pscoast -Dh  $bounds1 -O -K Lf-125/37.3/ -S100>> CaliforniaCalaverasMap1.ps



# start to plot the dem for the dry areas
# calcualte the gradient
grdgradient ttopo30.grd -V -A90  -Gtemp.grd
# illuminate the grid
grdmath temp.grd 9000. / = gtopo_grad.grd
# plot the land
grdimage $fig1 ttopo30.grd -Igtopo_grad.grd $bounds1 $proj1 -Cmydata.cpt -O -K -P >>  CaliforniaCalaverasMap1.ps
\rm temp.grd

# plot the coastline
gmtset LABEL_FONT_SIZE 12p   # makes the km label smaller
pscoast $fig1 $bounds1 $proj1  -W3 -A20 -Dh -O -P -K -Lf-122.15/39.75/39.75/50k >>  CaliforniaCalaverasMap1.ps
gmtset LABEL_FONT_SIZE 24p  # now you must remove it so the other labels are not so small


#switch off the clipping of the dry areas
pscoast $fig1 -Q -O -K $bounds1 $proj1 >>  CaliforniaCalaverasMap1.ps



# create the white filled box within which to place inset
#echo "236 36" > poly.txt
#echo "237.5 36" >> poly.txt
#echo "237.5 37.2" >> poly.txt
#echo "236 37.2" >> poly.txt
#psxy $fig1 $bounds1 $proj1 poly.txt -L -O -P -K -Gwhite>>  CaliforniaCalaverasMap1.ps


# join right side of insert to  area specified on large map
#echo "237.5     37.2" > poly2.txt
#echo "-121.6620  37.2847" >> poly2.txt
#echo "-121.6620  37.2945" >> poly2.txt
#echo "237.5     36" >> poly2.txt
#psxy $fig1 $bounds1 $proj1 poly2.txt -O -P -K  -Wblack >>  CaliforniaCalaverasMap1.ps


# join left side of insert to  area specified on large map
#echo "-121.6620  37.2945" >> poly4.txt
#echo "236 36" > poly4.txt
#echo "-121.6620  37.2945" >> poly4.txt
#echo "236 37" >> poly4.txt
#psxy $fig1 $bounds1 $proj1 poly4.txt -O -P -K  -Wwhite >>  CaliforniaCalaverasMap1.ps

# plot stations
awk '{ print $3 "  " $2}' station.dat   | psxy $fig1  $proj1 $bounds1   -St0.08i -O  -P -K -Gblack >>  CaliforniaCalaverasMap1.ps

# add a star at the location of the event cluster
echo "-121.6620  37.2847" >> poly3.txt
#psxy $fig1 $bounds1 $proj1 poly3.txt -Sa0.15i -O -P -K -Gblack >>  CaliforniaCalaverasMap1.ps
psxy $fig1 $bounds1 $proj1 poly3.txt -Sa0.2i -O -P -K  -Gwhite >>  CaliforniaCalaverasMap1.ps

gmtset ANNOT_FONT_SIZE 8p
# plot the events of interest with a red circle
#awk '{ print $5 "  " $6 "  " $17}' ../../hypoDD_loc68.txt | psxy $fig2 $anot2  $bounds2 $proj2 -Sc0.02i -O -P -Gred >>  CaliforniaCalaverasMap1.ps

gmtset ANNOT_FONT_SIZE 8p
# plot the events with a black circle
#awk '{ print $5 "  " $6 "  " $17}' hypoDD.reloc | psxy $fig2 $anot2  $bounds2 $proj2 -Sc0.02i -O -P -Gblack >>  CaliforniaCalaverasMap1.ps
#awk '{ print $5 "  " $6 "  " $17}' ../../hypodd_not68.txt | psxy $fig2 $anot2  $bounds2 $proj2 -Sc0.02i -P -Gblack >>  CaliforniaCalaverasMap1.ps


# Plot a polygon surrounding the events to be used. 
#echo "-577.4390    550.6190 " > poly4.txt
#echo " 31.0976    -477.7166 " >> poly4.txt
#echo "102.1098    -406.7044 " >> poly4.txt
#echo "-506.4268    621.6312 " >>poly4.txt
#echo "-577.4390    550.6190 " >> poly4.txt
#awk '{ print $1 "  " $2}' poly4.txt | psxy $fig2 $bounds2 $proj2 -P >>  CaliforniaCalaverasMap1.ps   #-L -O -P -K  >> tmp.ps

# plot the colourmap
#psscale $fig1 -Cdepth.cpt -D236/36/1/0.2h -O -P    >> CaliforniaCalaverasMap1.ps   
#-Ba10f5::/:M: -D15/3/6/0.6 -O  -P >> CaliforniaCalaverasMap1.ps

ps2epsi  CaliforniaCalaverasMap1.ps
eps2eps CaliforniaCalaverasMap1.epsi CaliforniaCalaverasMap1.eps
rm  CaliforniaCalaverasMap1.ps
rm CaliforniaCalaverasMap1.epsi
gv  CaliforniaCalaverasMap1.eps

rm .gmtdefaults4
