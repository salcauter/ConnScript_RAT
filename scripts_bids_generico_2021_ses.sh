#!/bin/bash

sourceDir=$1

############# CHECK for the parameters, files and PATH specified in this SECTION ###################

TR=1

anat=FLASH;
sigma=4.5;

#sub=sub-R020;
sub=*
#ses=ses-P91;
ses=*
run=*;

#sub=sub-R01_P120

# In Hz
hpf=0.01 
lpf=0.08

# Reference Volume for motion correction
refvol=0;

# ATLAS FILES  ### Check for the correct altasdir and atlas files !!! and their PATH
# atlasdir=/mnt/Data/RAT_fMRI/Templates_Rat/TohokuUniv/correctOrientation;
atlasdir=/misc/cannabis/alcauter/Desarrollo_RATs/TohokuUniv/correctOrientation/
atlasbrain=${atlasdir}/brain.nii.gz;
atlasbrain5mm=${atlasdir}/brain_5mm.nii.gz;
atlasbrain5mmMask=${atlasdir}/brain_5mm_mask.nii.gz;

# Create and populate the folder "derivatives"
mkdir -p ${sourceDir}/derivatives/ppBOLD
ls -d ${sourceDir}/${sub} | parallel cp -r {} ${sourceDir}/derivatives/ppBOLD


#  Voxel size (~x10)  ###  CORRECT FOR YOUR DATA !!!. All dimensions must be > 1, but keep the scale
# Corrige orientacion 
#echo ""
#echo "fslreorient2std anats"
#echo ""
#ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/*_${anat}.nii.gz | parallel fslreorient2std {}
#echo ""
#echo "fslreorient2std funcs"
#echo ""
#ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*gz | parallel fslreorient2std {}
echo ""
echo "fslreorient2std anats"
echo ""
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/*_${anat}.nii.gz | parallel fslorient -deleteorient {}
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/*_${anat}.nii.gz | parallel fslswapdim {} -x z -y {}
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/*_${anat}.nii.gz | parallel fslorient -setqformcode 1 {}
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/*_${anat}.nii.gz | parallel fslorient -setsformcode 1 {}
echo ""
echo "fslreorient2std funcs"
echo ""
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*rest*gz | parallel fslorient -deleteorient {}
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*rest*gz | parallel fslswapdim {} x -y z {}
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*rest*gz | parallel fslorient -setqformcode 1 {}
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*rest*gz | parallel fslorient -setsformcode 1 {}

echo ""
echo ">> Modifying VOXEL SIZE"
echo ""

for fl in $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*gz);do
	fx=$(fslval $fl pixdim1 | awk '{print $1*10}');
	fy=$(fslval $fl pixdim2 | awk '{print $1*10}');
	fz=$(fslval $fl pixdim3 | awk '{print $1*10}');

	fslchpixdim $fl $fx $fy $fz
	fslorient -setqformcode 1 $fl
	fslorient -setsformcode 1 $fl
done
for fl in $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/*_${anat}.nii.gz);do
	Tx=$(fslval $fl pixdim1 | awk '{print $1*10}');
	Ty=$(fslval $fl pixdim2 | awk '{print $1*10}');
	Tz=$(fslval $fl pixdim3 | awk '{print $1*10}');

	fslchpixdim $fl $Tx $Ty $Tz
	fslorient -setqformcode 1 $fl
	fslorient -setsformcode 1 $fl
done

############################  PREPROCESSING ########################################################

############# Slice timing and motion correction
echo ">> Slice timing and motion correction"
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/sub-*_task-rest_*.nii.gz | cut -d . -f 1 | parallel slicetimer -i {} -o {}_pp --odd 
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/sub-*_task-rest_*_pp.nii.gz | parallel mcflirt -in {} -refvol $refvol -plots

mkdir -p ${sourceDir}/derivatives/QC_pp
slicesdir $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/sub-*_task-rest_*_pp.nii.gz)
mv slicesdir ${sourceDir}/derivatives/QC_pp/slicesdir_func

############# CORREGISTROS    
echo ">> Structural preprocessing"
### BET Estructural
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}.nii.gz | cut -d . -f 1 | parallel N4BiasFieldCorrection -i {}.nii.gz -o {}_pp.nii.gz; # -t [0.3,0.01,100];
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}_pp.nii.gz | cut -d . -f 1 | parallel DenoiseImage -i {}.nii.gz -o {}.nii.gz;
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}.nii.gz | cut -d . -f 1 | parallel fslcpgeom {} {}_pp;
#ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}_pp.nii.gz | parallel fslswapdim {} x y z {}

for i in $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}_pp.nii.gz);do
	echo "Processing " $(basename $i .nii.gz | cut -d _ -f 1); 

	Tx=$(fslval $i pixdim1);
	Ty=$(fslval $i pixdim2);
	Tz=$(fslval $i pixdim3);
	Tzb=$(fslval $i pixdim3 | awk '{print $1/2}');


	fslchpixdim  $i $Tx $Ty $Tzb
	j=$(echo $i | cut -d . -f 1);
	center=$(cluster -i $i -t 10 | sed -n '2{p;q}' | awk '{print int($7)" "int($8)" "int($9)}');
	echo $i $j $center
	bet ${i} ${j}_brain.nii.gz -r 95 -f 0.25 -g 0.2; # -c $center ;
	bet ${j}_brain.nii.gz ${j}_brain.nii.gz -r 95 -f 0.25 -g -0.25; # -c $center 

	fslchpixdim $i $Tx $Ty $Tz
	fslchpixdim ${j}_brain.nii.gz $Tx $Ty $Tz

done

slicesdir -o $(ls -r ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}_pp*.nii.gz)
rm -rf ${sourceDir}/derivatives/QC_pp/slicesdir_BET; 
mv slicesdir/ ${sourceDir}/derivatives/QC_pp/slicesdir_BET

### anat brain to ATLAS
echo ">> Corregistration: anatbrain to atlas"
# register anatbrain to atlas
ls  ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}_pp_brain.nii.gz | cut -d . -f 1 | parallel flirt -in {} -ref $atlasbrain -out {}_2TohokuA -omat {}_2TohokuA.mat -bins 256 -cost corratio -searchrx -15 15 -searchry -15 15 -searchrz -15 15 -dof 12 -interp trilinear;

# Checar corregistros al atlas
slicesdir -p $atlasbrain $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}_pp_brain_2TohokuA.nii.gz)
rm -rf ${sourceDir}/derivatives/QC_pp/slicesdir_anattoAtlas; mv slicesdir/ ${sourceDir}/derivatives/QC_pp/slicesdir_anattoAtlas

### rsfMRI to anat
echo ">> Corregistration: rsfMRI to anat"
# reorient rsfMRIs
#ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/sub-*_task-rest_*_pp.nii.gz | cut -d . -f 1 | parallel fslswapdim {} x -y -z {}
# ExampleFunc
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/sub-*_task-rest_*_pp.nii.gz | cut -d . -f 1 | parallel fslroi {} {}_examplefunc $refvol 1;
# rsfMRI to anat
for i in $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/anat/sub-*_${anat}_pp_brain.nii.gz);do
	id=$(echo $i | rev | cut -d / -f 4 | rev);
	echo "Processing "$i;
	for j in $(ls ${sourceDir}/derivatives/ppBOLD/${id}/${ses}/func/${id}*_examplefunc.nii.gz | cut -d . -f 1);do 
	echo $(basename $j);

	N4BiasFieldCorrection -i ${j}.nii.gz -o ${j}NB.nii.gz;
	fslcpgeom ${j}.nii.gz ${j}NB.nii.gz;
	mv ${j}NB.nii.gz ${j}.nii.gz
	
	Fx=$(fslval $j pixdim1);
	Fy=$(fslval $j pixdim2);
	Fz=$(fslval $j pixdim3);
	Fzb=$(fslval $j pixdim3 | awk '{print $1/2}');


	fslchpixdim  $j $Fx $Fy $Fzb
	bj=$(echo $i | cut -d . -f 1);
	bet ${j} ${bj}_brain.nii.gz -r 95 -f 0.55 -g 0.2; # -c $center ;

	fslchpixdim $j $Fx $Fy $Fz
	fslchpixdim ${bj}_brain.nii.gz $Fx $Fy $Fz


	${FSLDIR}/bin/flirt -in ${bj}_brain -ref $i -out ${j}_2anat -omat ${j}_2anat.mat -cost corratio -dof 6 -interp trilinear -nosearch ;#-searchrx -15 15 -searchry -15 15 -searchrz -90 90 ;#-2D 
	done
done

# Slicesdir
rm ${sourceDir}/derivatives/listcheck_exfunc2anat.list
for i in $(ls -d ${sourceDir}/derivatives/ppBOLD/${sub});do 
	id=$(basename $i);
	for Ss in $(ls -d ${sourceDir}/derivatives/ppBOLD/${sub}/*);do
		echo ${sourceDir}/derivatives/ppBOLD/${id}/${Ss}/func/${id}*_task-rest_*_pp_examplefunc_2anat.nii.gz >> ${sourceDir}/derivatives/listcheck_exfunc2anat.list;
		echo ${sourceDir}/derivatives/ppBOLD/${id}/${Ss}/anat/${id}*_${anat}.nii.gz >> ${sourceDir}/derivatives/listcheck_exfunc2anat.list;
	done;
done
slicesdir -o $(cat ${sourceDir}/derivatives/listcheck_exfunc2anat.list)
rm -rf ${sourceDir}/derivatives/QC_pp/slicesdir_EF2anatw
mv slicesdir ${sourceDir}/derivatives/QC_pp/slicesdir_EF2anatw

### rsfMRI to ATLAS

# Combine transformations and apply
echo ">> Corregistration: rsfMRI to atlas"
for i in $(ls -d ${sourceDir}/derivatives/ppBOLD/${sub});do 
	echo $i;
	id=$(basename $i);	

	for pSs in $(ls -d ${sourceDir}/derivatives/ppBOLD/${id}/*);do
	Ss=$(basename $pSs);
	echo $Ss;

	anat_2TohA=$(ls ${sourceDir}/derivatives/ppBOLD/${id}/${Ss}/anat/${id}_*_${anat}_pp_brain_2TohokuA.mat);
	for j in $(ls ${sourceDir}/derivatives/ppBOLD/${id}/${Ss}/func/${id}*_EPI.nii.gz | cut -d . -f 1);do 
		#combine
		${FSLDIR}/bin/convert_xfm -omat ${j}_pp_examplefunc_2atlasbrain.mat $anat_2TohA -concat  ${j}_pp_examplefunc_2anat.mat;
		# Apply
		${FSLDIR}/bin/flirt -in ${j}_pp_mcf.nii.gz -applyxfm -init ${j}_pp_examplefunc_2atlasbrain.mat -out ${j}_pp_mcf_2atlasbrain5mm -paddingsize 0.0 -interp trilinear -ref $atlasbrain5mm;
		done	
	done
done

slicesdir -p $atlasbrain5mm $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/${Ss}/func/*_pp_mcf_2atlasbrain5mm.nii.gz)
rm -rf ${sourceDir}/derivatives/QC_pp/slicesdir_ppAtlas5mm
mv slicesdir ${sourceDir}/derivatives/QC_pp/slicesdir_ppAtlas5mm

############# aCompCor and Band Passing 
echo ">> aCompCor and Band Passing";
for i in $(ls -d ${sourceDir}/derivatives/ppBOLD/${sub});do 
	echo $i;
	id=$(basename $i);

	for pSs in $(ls -d ${sourceDir}/derivatives/ppBOLD/${id}/${ses});do
	Ss=$(basename $pSs);
	echo $Ss;

		for j in $(ls ${sourceDir}/derivatives/ppBOLD/${id}/${Ss}/func/${id}*_EPI_pp_mcf_2atlasbrain5mm.nii.gz | cut -d . -f 1);do
			liminf=$(fslstats $j -r | awk '{print $1}');
			liminf2=$(echo "scale=2; 2*${liminf}" | bc);
			${FSLDIR}/bin/fslmaths $j -Tmean -thr $liminf2 -bin -mas ${atlasdir}/CSFWMmask_5mm.nii.gz ${j}_CSFWMmask_5mm.nii.gz
			${FSLDIR}/bin/fslmeants -i $j -o ${j}_aCompCor.txt -m ${j}_CSFWMmask_5mm.nii.gz --eig --order=5;

			rm -f ${j}_aCCtf.nii.gz
			3dTproject -input ${j}.nii.gz -prefix ${j}_aCCtf.nii.gz -polort 0 -ort ${j}_aCompCor.txt -passband $hpf $lpf -TR ${TR} -mask $atlasbrain5mmMask;
		done
	done
done

### Smoothing 
echo ">> Smoothing"
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*EPI_pp_mcf_2atlasbrain5mm_aCCtf.nii.gz | cut -d . -f 1 | parallel fslmaths {} -s $sigma -mas $atlasbrain5mmMask {}_sm.nii.gz

echo ">> Connectivity Maps: ACC"
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*EPI_pp_mcf_2atlasbrain5mm_aCCtf_sm.nii.gz | cut -d . -f 1 | parallel fslmeants -i {} -o {}_RSCx.ts -m ${atlasdir}/brain_5mm_rscx.nii.gz
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*EPI_pp_mcf_2atlasbrain5mm_aCCtf_sm.nii.gz | cut -d . -f 1 | parallel 3dTcorr1D -prefix {}_RSCx.nii.gz -mask ${atlasdir}/brain_5mm_mask.nii.gz {}.nii.gz {}_RSCx.ts

echo ">> Connectivity Maps: M1"
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*EPI_pp_mcf_2atlasbrain5mm_aCCtf_sm.nii.gz | cut -d . -f 1 | parallel fslmeants -i {} -o {}_M1.ts -m ${atlasdir}/M1_mask_ero_5mm.nii.gz
ls ${sourceDir}/derivatives/ppBOLD/${sub}/${ses}/func/*EPI_pp_mcf_2atlasbrain5mm_aCCtf_sm.nii.gz | cut -d . -f 1 | parallel 3dTcorr1D -prefix {}_M1.nii.gz -mask ${atlasdir}/brain_5mm_mask.nii.gz {}.nii.gz {}_M1.ts



############# MENSAJE FINAL #############

# firefox --new-window ${sourceDir}/derivatives/QC_pp/slicesdir_func/index.html;
# firefox --new-tab ${sourceDir}/derivatives/QC_pp/slicesdir_BET/index.html;
# firefox --new-tab ${sourceDir}/derivatives/QC_pp/slicesdir_anattoAtlas/index.html;
# firefox --new-tab ${sourceDir}/derivatives/QC_pp/slicesdir_EF2anatw/index.html;
# firefox --new-tab ${sourceDir}/derivatives/QC_pp/slicesdir_ppAtlas5mm/index.html;
echo ">>         DONE"
echo ""


############# Optional but RECOMMENDED to save space:
# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_*.nii.gz
# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_*_pp.nii.gz
# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_*_pp_mcf.nii.gz
# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_*_pp_mcf_2atlasbrain5mm.nii.gz



################## TimeFiltering and smoothing Alternatives (using FSL only) #########################

### Smoothing: DenoiseImage (takes about 2-3 hours per 4D series)
#ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/*EPI_pp_mcf_2atlasbrain5mm_aCCtf.nii.gz | cut -d . -f 1 | parallel DenoiseImage -i {} -x $atlasbrain5mmMask -o {}_smDI.nii.gz

### aCompCor and time filtering with FSL

# fsl_glm
# ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/*EPI_pp_mcf_2atlasbrain5mm.nii.gz | cut -d . -f 1 | parallel fsl_glm -i {} -m $atlasdirbrain5mmMask -d {}_aCompCor.txt --demean --out_res={}_aCC.nii.gz;

# TR=1;
# hpf=0.01;
# lpf=0.1;
# hp_sigma=`echo "scale=2 ;(1/${hpf})/2.35/${TR}" | bc`; # In volumes for fslmaths
# lp_sigma=`echo "scale=2 ;(1/${lpf})/2.35/${TR}" | bc`; # In volumes for fslmaths

# ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/*EPI_pp_mcf_2atlasbrain5mm_aCC.nii.gz | cut -d . -f 1 | parallel fslmaths {} -bptf  $hp_sigma $lp_sigma -mas $atlasbrain5mmMask {}tf.nii.gz

# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/*EPI_pp_mcf_2atlasbrain5mm_aCC.nii.gz

