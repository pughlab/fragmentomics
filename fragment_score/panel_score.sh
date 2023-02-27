INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_TS/all_unique
PDIR=/cluster/projects/pughlab/projects/CHARM/LFS/fragment_score
shdir=$PDIR/sh_scripts/panel
outdir=$PDIR/output/panel
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
module load samtools\n
module load R/4.0.0\n" > $shdir/${bam}.sh

echo -e "samtools view -bh $INPUTDIR/${bam}.bam chr17:7668000-7688000 > $outdir/${bam}.bam
samtools index $outdir/${bam}.bam\n" >> $shdir/${bam}.sh

echo -e "Rscript $scripts/04_patient_score.R\
 --id $bam\
 --bam $outdir/${bam}.bam\
 --ref $ref\
 --outdir $outdir\
 --libdir $scripts\n" >> $shdir/${bam}.sh

echo -e "if [[ -f \"$outdir/${bam}*score.txt\" ]]
then
rm $outdir/${bam}.bam*" >> $shdir/${bam}.sh

done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 24:00:00 -p all --mem 16G $file
done
