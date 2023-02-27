INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_TS/all_unique
PDIR=/cluster/projects/pughlab/projects/CHARM/LFS/fragment_score
germline=/cluster/projects/pughlab/projects/CHARM/LFS/pipeline_output/HaplotypeCaller/cohort/germline_variants
somatic=/cluster/projects/pughlab/projects/CHARM/LFS/pipeline_output/MuTect2
oicr=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_TS/vcfs
shdir=$PDIR/sh_scripts/variant
outdir=$PDIR/output/variant
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
name=${bam:0:10}
id=$(echo $bam | sed 's/\_TS.*/_TS/')
echo $bam
echo $name
echo $id

g_vcf=$(ls -t $germline/${name}*vcf.gz)
s_vcf=$(ls -t $somatic/$name/${id}_all/*vcf.gz)
o_vcf=$(ls -t $oicr/${id}*all.unique*vcf.gz)
echo $g_vcf
echo $s_vcf
echo $o_vcf

echo -e "#!/bin/bash\n
module load R/4.0.0\n" > $shdir/${bam}.sh

echo -e "Rscript $scripts/05_variant_score.R\
 --id $id\
 --bam $INPUTDIR/${bam}.bam\
 --ref $ref\
 --vcf $g_vcf\
 --type germline\
 --outdir $outdir\
 --libdir $scripts\n" >> $shdir/${bam}.sh

echo -e "Rscript $scripts/05_variant_score.R\
 --id $id\
 --bam $INPUTDIR/${bam}.bam\
 --ref $ref\
 --vcf $s_vcf\
 --type somatic\
 --outdir $outdir\
 --libdir $scripts" >> $shdir/${bam}.sh

echo -e "Rscript $scripts/05_variant_score.R\
 --id $id\
 --bam $INPUTDIR/${bam}.bam\
 --ref $ref\
 --vcf $o_vcf\
 --type oicr\
 --outdir $outdir\
 --libdir $scripts" >> $shdir/${bam}.sh

done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 24:00:00 -p all --mem 8G $file
done
