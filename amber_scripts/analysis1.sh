module load amber/20_ompigw
export OMP_NUM_THREADS=1
ulimit -s unlimited

NAME=$(head -n1 sub_history | cut -d "/" -f 2 | cut -d " " -f 1)
FOL=$(pwd | cut -d "_" -f 2)
mkdir -p analysis_$FOL && cd analysis_$FOL

if [[ $NAME == *"P50416"* || $NAME == *"2H3WA"* || $NAME == *"CPT1A"* ]]; then
   RESI=(82 314 530 553)
   RESN=(TYR241 HID473 THR689 PHE712)
else
   RESI=(90 307 523 546)
   RESN=(TYR256 HID473 ILE689 PHE712)
fi

#
# MIN
#
#process_mdout.perl ../min*.out
#mv summary.TEMP min.TEMP
#mv summary.PRES min.PRES
#mv summary.EKTOT min.EKTOT
#mv summary.EPTOT min.EPTOT
#mv summary.ETOT min.ETOT
#mv summary.DENSITY min.DENS
#mv summary.VOLUME min.VOL
#rm summary*

#
#EQUIL and PROD
#mpirun -np 6 cpptraj.MPI <<EOF
#   readdata ../prod*.out name MDOUT
#   writedata prod.ETOT  MDOUT[Etot]    time 0.0002
#   writedata prod.EPTOT MDOUT[EPtot]   time 0.0002
#   writedata prod.EKTOT MDOUT[EKtot]   time 0.0002
#   writedata prod.VDW   MDOUT[VDW]     time 0.0002
#   writedata prod.EEL   MDOUT[EEL1-4]  time 0.0002
#   writedata prod.DENS  MDOUT[Density] time 0.0002
#   writedata prod.TEMP  MDOUT[TEMP]    time 0.0002
#   writedata prod.VOL   MDOUT[VOLUME]  time 0.0002
#   writedata prod.PRES  MDOUT[PRESS]   time 0.0002
#   run
#   clear all
#EOF

#mpirun -np 1 cpptraj.MPI <<EOF
#   parm ../*parm7
#   trajin ../prod*.mdcrd 1 last 100
#   autoimage
#   strip :WAT,Na+,Cl-
#   rms first mass @C,CA,N
#   trajout traj_${NAME}_FULL.pdb offset 2
#   trajout traj_${NAME}_0500.pdb onlyframes 1-500
#   trajout traj_${NAME}_1000.pdb onlyframes 500-1000
#   trajout traj_${NAME}_1500.pdb onlyframes 1000-1500
#   trajout traj_${NAME}_2000.pdb onlyframes 1500-2000
#   trajout traj_${NAME}_2500.pdb onlyframes 2000-2500
#   trajout traj_${NAME}_3000.pdb onlyframes 2500-3000
#   run
#   clear all
#EOF

# CALCULATE RMSF
mpirun -np 6 cpptraj.MPI <<EOF
# CREATE STRIPPED TOPOLOGY
#   parm ../com_solv.parm7
#   parmstrip :WAT,Na+,Cl-
#   parmwrite out ../com_gas.parm7
#   run
#   clear all

# CREATE FULL TRAJ
#   parm ../*parm7
#   trajin ../prod*ns.mdcrd 1 last 10
#   autoimage
#   rms first mass @C,CA,N
#   strip :WAT,Na+,Cl-
#   trajout ../prod_FULL.mdcrd
#   run
#   clear all

# CREATE AVG STRUCTURE
#   parm ../com_gas.parm7 
#   trajin ../prod_FULL.mdcrd 14000 last
#   autoimage
#   rms first mass @C,CA,N
#   average avg_100ns.rst7
#   average avg_100ns.pdb
#   rms ref MyAvg
#   atomicfluct out backbone.RMSF @C,CA,N&!@H= byres
#   run

   parm ../com_solv.parm7
   reference ../prod1500ns.rst7 lastframe [MyRef] 
   trajin ../prod*.mdcrd 1 last 200
   strip :WAT,Na+,Cl-
#   average crdset MyAvg start 2500 stop 3000

  rms Backbone1 ref [MyRef] :10-119,136-146,167-181,216-254,271-323,361-500,518-525,541-590@C,N,CA time 2 out backbone_noloops.RMSD
 #  rms Backbone2 @C,N,O,CA&!@H= ref MyAvg time 2 out backbone_all.RMSD
 #  rms Pocket :MOL<:4.0&!(:MOL) ref [MyRef] mass out pocket.RMSD
 #  rms Ligand :MOL nofit ref [MyRef] mass out ligand.RMSD

 #  multidihedral Chi1 chip chi2 resrange ${RESI[0]}-${RESI[0]} out ${RESN[0]}.dat
 #  multidihedral Chi1 chip chi2 resrange ${RESI[1]}-${RESI[1]} out ${RESN[1]}.dat
 #  multidihedral Chi1 chip chi2 resrange ${RESI[2]}-${RESI[2]} out ${RESN[2]}.dat
 #  multidihedral Chi1 chip chi2 resrange ${RESI[3]}-${RESI[3]} out ${RESN[3]}.dat

#  distance DIST :95@OH :313@O out TYR254-dist1.dat
   run
EOF

