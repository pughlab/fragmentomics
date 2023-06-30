#!/bin/bash

analysis=hematopoetic
CPU=1
mem=8G

griffin=/cluster/projects/pughlab/bin/Griffin/v0.2.0
basedir=/cluster/projects/pughlab/projects/CHARM/LFS/griffin2
sites=$griffin/site_configs/${analysis}_sites.yaml
ref=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
input=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
counts=$basedir/output/GC_correction
outdir=$basedir/output/nucleosome_profiling/$analysis
shdir=$basedir/sh_scripts/nucleosome_profiling/$analysis

encode_exclude=$griffin/Ref/encode_unified_GRCh38_exclusion_list.bed
centromere_path=$griffin/Ref/hg38_centromeres.bed
gap_path=$griffin/Ref/hg38_gaps.bed
patch_path=$griffin/Ref/hg38_fix_patches.bed
alternative_haplotype_path=$griffin/Ref/hg38_alternative_haplotypes.bed

mkdir -p $outdir
mkdir -p $outdir/tmp
mkdir -p $outdir/results
mkdir -p $shdir

cd $input
ls *bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
mv bam bams

for bam in $(cat bams);do

name=${bam:0:25}
echo $bam
echo $name

echo -e "#!/bin/bash\n
source activate base\n
conda activate griffin2\n" > $shdir/${name}.sh

echo -e "$griffin/scripts/griffin_coverage.py \
--sample_name $name \
--bam $input/${bam}.bam \
--GC_bias $counts/GC_bias/${name}.GC_bias.txt \
--mappability_bias $counts/mappability_bias/${name}.mappability_bias.txt \
--mappability_correction True \
--tmp_dir $outdir/tmp \
--reference_genome $ref \
--mappability_bw $griffin/Ref/k50.Umap.MultiTrackMappability.hg38.bw \
--chrom_sizes_path $griffin/Ref/hg38.standard.chrom.sizes \
--sites_yaml $sites \
--griffin_scripts $griffin/scripts \
--chrom_column Chrom \
--position_column position \
--strand_column Strand \
--chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
--norm_window -5000 5000 \
--size_range 100 200 \
--map_quality 20 \
--number_of_sites none \
--sort_by none \
--ascending none \
--CPU $CPU\n" >> $shdir/${name}.sh

echo -e "$griffin/scripts/griffin_merge_sites.py \
--sample_name $name \
--uncorrected_bw_path $outdir/tmp/$name/tmp_bigWig/${name}.uncorrected.bw \
--GC_corrected_bw_path $outdir/tmp/$name/tmp_bigWig/${name}.GC_corrected.bw \
--GC_map_corrected_bw_path $outdir/tmp/$name/tmp_bigWig/${name}.GC_map_corrected.bw \
--mappability_correction False \
--tmp_dir $outdir/tmp \
--results_dir $outdir/results \
--mappability_bw $griffin/Ref/k50.Umap.MultiTrackMappability.hg38.bw \
--chrom_sizes_path $griffin/Ref/hg38.standard.chrom.sizes \
--sites_yaml $sites \
--griffin_scripts $griffin/scripts \
--chrom_column Chrom \
--position_column position \
--strand_column Strand \
--chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
--norm_window -5000 5000 \
--save_window -1000 1000 \
--fft_window -960 960 \
--fft_index 10 \
--smoothing_length 165 \
--exclude_paths $encode_exclude $centromere_path $gap_path $patch_path $alternative_haplotype_path \
--step 15 \
--CNA_normalization False \
--individual False \
--smoothing True \
--exclude_outliers True \
--exclude_zero_mappability True \
--number_of_sites none \
--sort_by none \
--ascending none \
--CPU $CPU\n" >> $shdir/${name}.sh

done

cd $shdir
ls *.sh > files
for file in $(cat files);do
sbatch -c $CPU --mem $mem -t 24:00:00 $file
done
