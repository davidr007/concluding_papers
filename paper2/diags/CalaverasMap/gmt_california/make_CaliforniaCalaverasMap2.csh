#!/bin/csh

gmtset PAPER_MEDIA A4

#set fig1 = "-Xa1i -Ya1i"
#set bounds1 = "-R-124/-120/36/40"
#set proj1 = "-JM5i"
#set anot1 = " -B1g30"




#set bounds2 = "-R-121.8/-121.5/37.1/37.4"
set bounds2 = "-R-600/600/-600/800"
set proj2 = "-JX10c/11.167c"
#set anot2 = " -BnSeW400"
set fig2 = "-Xa1.6i -Ya1.5i"

# COLOR_MODEL = RGB
echo "0 215 025 028  1 215 025 028" > depth.cpt
echo "1 253 174 097  2 253 174 097" >> depth.cpt
echo "2 255 255 191  3 255 255 191" >> depth.cpt
echo "3 171 221 164  4 171 221 164" >> depth.cpt
echo "4 043 131 186  5 043 131 186 " >> depth.cpt
# plot the events by magnitude
awk '{ print $5 "  " $6 "  " $17}' hypoDD.reloc | psxy $fig2 -BnSeWa400f100:"x (m)":/a400f100:"y (m)":  $bounds2 $proj2 -Sc0.04i -K -P -Cdepth.cpt >!  CaliforniaCalaverasMap2.ps

# Plot a polygon surrounding the events to be used. 
echo "-577.4390    550.6190 " > poly3.txt
echo " 31.0976    -477.7166 " >> poly3.txt
echo "102.1098    -406.7044 " >> poly3.txt
echo "-506.4268    621.6312 " >>poly3.txt
echo "-577.4390    550.6190 " >> poly3.txt
awk '{ print $1 "  " $2}' poly3.txt | psxy $fig2 $bounds2 $proj2 -P -O -K -W0.04c>>  CaliforniaCalaverasMap2.ps   #-L -O -P -K  >> tmp.ps

# plot the colourmap
#psscale $fig2 -Cdepth.cpt -D236/36/1/0.2h -O -P    >> CaliforniaCalaverasMap2.ps   
psscale $fig2 -D8c/8c/4c/0.6c -B1::/:"km": -O -P -Cdepth.cpt >>  CaliforniaCalaverasMap2.ps

ps2epsi  CaliforniaCalaverasMap2.ps
mv CaliforniaCalaverasMap2.epsi CaliforniaCalaverasMap2.eps
rm  CaliforniaCalaverasMap2.ps
rm CaliforniaCalaverasMap2.epsi
gv  CaliforniaCalaverasMap2.eps

