
out=log.lammps     # Archivo de salida de lammps
# Genera el archivo -- md.dat -- que pega toda la info termodinámica
# correspond. a cada rampa y meseta de temperatura. 

cat $out|sed -n "/Step/,/Loop/p"|head -n-1|tail -n+2  > DATA
sed -i '$ d' DATA		# erase last line in file

sed -i '/Loop/d' DATA
sed -i '/Step/d' DATA

head=$(grep "Step " $out| head -1)
echo "# $head"        > md.dat
cat DATA            >> md.dat 
rm -f DATA



# Genera el archivo -- EvsT -- listo para plotear en xmgrace
cat EvsT-lammps.dat | awk {'print $2"  " $3'} > EvsT.dat
sed -i '/Time/d' EvsT.dat
