NAME=$(head -n1 sub_history | cut -d "/" -f 2 | cut -d " " -f 1)
mkdir -p analysis_$NAME && cd analysis_$NAME

if [[ $NAME == *"P50416"* ]] || [[ $NAME == *"2H3WA"* ]]; then
   RESI=(82 314 530 553)
   RESN=(TYR241 HID473 THR689 PHE712)
else
   RESI=(90 307 523 546)
   RESN=(TYR256 HID473 ILE689 PHE712)
fi

#
# MIN
#
process_mdout.perl ../min*.out
mv summary.TEMP min.TEMP
mv summary.PRES min.PRES
mv summary.EKTOT min.EKTOT
mv summary.EPTOT min.EPTOT
mv summary.ETOT min.ETOT
mv summary.DENSITY min.DENS
mv summary.VOLUME min.VOL
rm summary*

#
#EQUIL and PROD
#
cpptraj <<EOF
   readdata ../prod*.out name MDOUT
   writedata prod.ETOT  MDOUT[Etot]    time 0.0002
   writedata prod.EPTOT MDOUT[EPtot]   time 0.0002
   writedata prod.EKTOT MDOUT[EKtot]   time 0.0002
   writedata prod.VDW   MDOUT[VDW]     time 0.0002
   writedata prod.EEL   MDOUT[EEL1-4]  time 0.0002
   writedata prod.DENS  MDOUT[Density] time 0.0002
   writedata prod.TEMP  MDOUT[TEMP]    time 0.0002
   writedata prod.VOL   MDOUT[VOLUME]  time 0.0002
   writedata prod.PRES  MDOUT[PRESS]   time 0.0002
   run
EOF

cpptraj <<EOF

   parm ../*parm7
   reference ../prod*10ns.rst7 [Last10] 
   for TRAJ in ../prod*.mdcrd
      trajin \$TRAJ 1 last 100
   done

   rms Backbone :10-595@C,N,O,CA ref [Last10] mass out backbone.RMSD
   rms Pocket :COA<:4.0&!(:WAT,Na+,Cl-,MOL) ref [Last10] mass out pocket.RMSD
   rms Ligand :COA nofit ref [Last10] mass out ligand.RMSD

   #for l=0;l<4;l++
   multidihedral Chi1 chip chi2 resrange ${RESI[0]}-${RESI[0]} out ${RESN[0]}.dat
   multidihedral Chi1 chip chi2 resrange ${RESI[1]}-${RESI[1]} out ${RESN[1]}.dat
   multidihedral Chi1 chip chi2 resrange ${RESI[2]}-${RESI[2]} out ${RESN[2]}.dat
   multidihedral Chi1 chip chi2 resrange ${RESI[3]}-${RESI[3]} out ${RESN[3]}.dat
   #done

#   distance DIST :95@OH :313@O out TYR254-dist1.dat

   run
EOF

