INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
PDIR=/cluster/projects/pughlab/projects/CHARM/LFS/fragment_score
shdir=$PDIR/sh_scripts/patient
outdir=$PDIR/output/patient
scripts=/cluster/projects/pughlab/bin/fragmentomics/v2/fragment_score/R
ref=/cluster/projects/pughlab/bin/fragmentomics/v2/fragment_score/ref/vessies_reference_set.txt

mkdir -p $shdir
mkdir -p $outdir

cd $INPUTDIR
ls *.bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
echo -e "#!/bin/bash\n
module load R/4.0.0\n" > $shdir/${bam}.sh

echo -e "Rscript $scripts/04_patient_score.R\
 --id $bam\
 --bam $INPUTDIR/${bam}.bam\
 --ref $ref\
 --outdir $outdir\
 --libdir $scripts" >> $shdir/${bam}.sh

done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 24:00:00 -p all --mem 16G $file
done
