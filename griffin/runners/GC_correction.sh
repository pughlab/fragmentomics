#!/bin/bash

griffin=/cluster/projects/pughlab/bin/Griffin/v0.2.0
basedir=/cluster/projects/pughlab/projects/CHARM/LFS/griffin2
ref=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
input=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
outdir=$basedir/output/GC_correction
shdir=$basedir/sh_scripts/GC_correction

mkdir -p $outdir
mkdir -p $shdir
mkdir -p $outdir/mappability_bias
mkdir -p $outdir/mappability_plots
mkdir -p $outdir/tmp

cd $input
ls *bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
mv bam bams

for bam in $(cat bams);do

name=${bam:0:25}
echo $bam
echo $name

mkdir -p $outdir/tmp/$name

echo -e "#!/bin/bash
source activate base
conda activate griffin2" > $shdir/${name}.sh

echo -e "$griffin/scripts/griffin_mappability_correction.py \
--bam_file $input/${bam}.bam \
--bam_file_name $name \
--output $outdir/mappability_bias/${name}.mappability_bias.txt \
--output_plot $outdir/mappability_plots/${name}.mappability_bias.pdf \
--mappability $griffin/Ref/k50.Umap.MultiTrackMappability.hg38.bw \
--exclude_paths $griffin/Ref/encode_unified_GRCh38_exclusion_list.bed \
--chrom_sizes $griffin/Ref/hg38.standard.chrom.sizes \
--map_quality 20 \
--CPU 8 \
--tmp_dir $outdir/tmp/$name" >> $shdir/${name}.sh

echo -e "$griffin/scripts/griffin_GC_counts.py \
--bam_file $input/${bam}.bam \
--bam_file_name $name \
--mappable_regions_path $griffin/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed \
--ref_seq $ref \
--chrom_sizes $griffin/Ref/hg38.standard.chrom.sizes \
--out_dir $outdir \
--map_q 20 \
--size_range 15 500 \
--CPU 8" >> $shdir/${name}.sh

echo -e "$griffin/scripts/griffin_GC_bias.py \
--bam_file_name $name \
--mappable_name k100_minus_exclusion_lists.mappable_regions.hg38 \
--genome_GC_frequency $griffin/Ref/genome_GC_frequency \
--out_dir $outdir/ \
--size_range 15 500" >> $shdir/${name}.sh

done

cd $shdir
ls *.sh > files
for file in $(cat files);do
sbatch -p all -c 8 --mem 8G -t 24:00:00 $file
done
