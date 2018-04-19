#!/bin/bash
# SCOPE: GraphCut with free parameters rej, fmin, fmax, wm. Normalized images. It runs on Igrida cluster.

T1masked=$1
T1Seg=$2
DWITensor=$3
Mask=$4
DWITracto=$5

echo $1
echo $4

export PATH=$PATH:~/Anima-Build/bin/

. /etc/profile.d/modules.sh

module load neurinfo/fsl
source $FSLDIR/etc/fslconf/fsl.sh

~/linux_openmp_64/3dresample -master $T1masked -prefix aparc+aseg_T1.nii.gz -inset $T1Seg -overwrite
if [ ! -f $DWITracto ]; then
    animaDTITractography --max-length 150 --min-length 10 --nb-fibers 1 -o $DWITracto -s $Mask -i $DWITensor -T 10
fi

COUNTER=1
while (( $COUNTER < 81 ))
do
	echo $COUNTER
	COUNTER2=$COUNTER
	let COUNTER2+=1
	mkdir Freesurfer_tract
	fslmaths aparc+aseg_T1.nii.gz -thr ${v[$COUNTER]} -uthr ${v[$COUNTER]} tmp_mask.nii.gz
	fslmaths tmp_mask.nii.gz -bin tmp_mask_bin.nii.gz

	animaFibersFilterer -t 1 -r tmp_mask_bin.nii.gz -i $DWITracto -o tmp_dwi_fiber_ROI.vtk -T 10

	mv tract_length.txt Freesurfer_tract
	mv tract_number.txt Freesurfer_tract

	mv Freesurfer_tract/tract_length.txt Freesurfer_tract/tract_length_${COUNTER}.txt
	mv Freesurfer_tract/tract_number.txt Freesurfer_tract/tract_number_${COUNTER}.txt

	while (( $COUNTER2 < 81 ))
	do
		echo $COUNTER2
		fslmaths aparc+aseg_T1.nii.gz -thr ${v[$COUNTER2]} -uthr ${v[$COUNTER2]} tmp_mask_2.nii.gz
		fslmaths tmp_mask_2.nii.gz -bin tmp_mask_2_bin.nii.gz

		fslmaths tmp_mask_2_bin.nii.gz -mul 2 tmp_mask22.nii.gz
		fslmaths tmp_mask_bin.nii.gz -add tmp_mask22.nii.gz tmp_mask_gobal.nii.gz

		animaFibersFilterer -t 1 -t 2 -r tmp_mask_gobal.nii.gz -i $DWITracto -o tmp_dwi_fiber_ROI.vtk _T 10

		mv tract_length.txt Freesurfer_tract
		mv tract_number.txt Freesurfer_tract

		mv Freesurfer_tract/tract_length.txt Freesurfer_tract/tract_length_${COUNTER}_${COUNTER2}.txt
		mv Freesurfer_tract/tract_number.txt Freesurfer_tract/tract_number_${COUNTER}_${COUNTER2}.txt

		let COUNTER2+=1
		echo $COUNTER2
	done
	let COUNTER+=1

done
rm -rf tmp*
echo "fin d'un patient"
echo $COUNTER2
echo $COUNTER
