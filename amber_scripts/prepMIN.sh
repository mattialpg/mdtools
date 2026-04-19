
# Look for C-terminus and water molecules
A="$(grep 'MOL' $1.pdb | tail -1)"
B="$(grep 'WAT' $1.pdb | head -1)"
C="$(grep 'WAT' $1.pdb | tail -1)"

printf "\n$A\n\n$B\n$C\n\n"

# Strip the lines
read -ra An -d ''<<<"$A"
read -ra Bn -d ''<<<"$B"
read -ra Cn -d ''<<<"$C"

# Replace and create new files
sed "s/resA/${An[4]}/g" ~/.bin/min1.in > min1.x
sed "s/resB/${Bn[4]}/g" ~/.bin/min2.in > min2.x
sed -i "s/resC/${Cn[4]}/g" min2.x
sed "s/resA/${An[4]}/g" ~/.bin/min3.in > min3.x
sed "s/resA/${An[4]}/g" ~/.bin/min4.in > min4.x
sed "s/resA/${An[4]}/g" ~/.bin/min5.in > min5.x
sed "s/resA/${An[4]}/g" ~/.bin/min6.in > min6.x
