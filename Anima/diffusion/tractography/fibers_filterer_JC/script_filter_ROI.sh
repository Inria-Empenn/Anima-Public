#!/bin/bash
# SCOPE: GraphCut with free parameters rej, fmin, fmax, wm. Normalized images. It runs on Igrida cluster.

# DATA dir
export PATH=$PATH:~/Anima-Build/bin/
#export PATH=$PATH:~/Dev/scripts/tracto_freesurfer/
# DATA dir
DATADIR=/temp_dd/igrida-fs1/jcoloign/longidep/Seg_freesurfer/
echo "$DATADIR"
# OUTPUT dir
OutputDirBase=/temp_dd/igrida-fs1/jcoloign/longidep/Seg_freesurfer/
mkdir -p $OutputDirBase
echo "$OutputDirBase"
#v=('1001' '1002' '1003' '1005' '1006' '1007' '1008' '1009' '1010' '1011' '1012' '1013' '1014' '1015' '1016' '1017' '1018' '1019' '1020' '1021' '1022' '1023' '1024' '1025' '1026' '1027' '1028' '1029' '1030' '1031' '1032' '1033' '1034' '1035' '2001' '2002' '2003' '2005' '2006' '2007' '2008' '2009' '2010' '2011' '2012' '2013' '2014' '2015' '2016' '2017' '2018' '2019' '2020' '2021' '2022' '2023' '2024' '2025' '2026' '2027' '2028' '2029' '2030' '2031' '2032' '2033' '2034' '2035' '10' '11' '12' '13' '17' '18' '49' '50' '51' '52' '53' '54')
#Nv=81
for patient in $DATADIR/*; do
   patientID=$(basename "$patient")
	 echo "$patientID"
	 pIn="$patient"
	 echo "$patient"
   cd $pIn
   patientID=$(basename "$patient")
   mkdir Freesurfer_tract
   rmdir Freesurfer_tract_
   chmod +x ~/Dev/scripts/tracto_freesurfer/main_tracto_ROI.sh
   #oarsub -l nodes=10,walltime=20:00:00 -p "cluster='neurinfo1'" "bash ~/Dev/scripts/tracto_freesurfer/main_tracto_ROI.sh /temp_dd/igrida-fs1/jcoloign/longidep/DTI_atlas/T13D_${patientID}_masked.nii.gz /temp_dd/igrida-fs1/jcoloign/longidep/Seg_freesurfer/${patientID}/aparc+aseg.nii.gz /temp_dd/igrida-fs1/jcoloign/longidep/DTI_atlas/dwi_${patientID}_Tensors.nii.gz /temp_dd/igrida-fs1/jcoloign/longidep/DTI_atlas/T13D_${patientID}_brainMask.nii.gz /temp_dd/igrida-fs1/jcoloign/longidep/Seg_freesurfer/${patientID}/dwi_${patientID}_fiber_Tensor.vtk"
done
