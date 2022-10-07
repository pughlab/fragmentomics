INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
basedir=/cluster/projects/pughlab/projects/CHARM/LFS/end_motifs
ref=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
picard_dir=/cluster/tools/software/picard/2.10.9
frags=/cluster/projects/pughlab/bin/fragmentomics/R

site=TP53

bed=/cluster/projects/pughlab/projects/CHARM/LFS/dinucleotide/${site}_regions.bed
outdir=$basedir/output/$site
shdir=$basedir/sh_scripts/$site

mkdir -p $shdir
mkdir -p $outdir

cd $INPUTDIR
ls *.bam > $shdir/bams
cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
### Load in the modules
echo -e "#!/bin/bash\n
module load samtools/1.10
module load picard/2.10.9
module load bedtools/2.27.1
module load R/4.0.0\n" > $shdir/${bam}.sh

### Make bam file of only site regions
echo -e "samtools view -hb -f 0x2 -L $bed $INPUTDIR/${bam}.bam > $outdir/${bam}_site.bam
samtools index $outdir/${bam}_site.bam\n" >> $shdir/${bam}.sh

### Remove duplicates
echo -e "java -jar $picard_dir/picard.jar MarkDuplicates \
I=$outdir/${bam}_site.bam \
O=$outdir/${bam}_deduped.bam \
M=$outdir/${bam}_metrics.txt \
TMP_DIR=$outdir \
REMOVE_DUPLICATES=true \
REMOVE_SEQUENCING_DUPLICATES=true

samtools sort -n $outdir/${bam}_deduped.bam -o $outdir/${bam}_deduped_sorted.bam\n" >> $shdir/${bam}.sh

### Sort bam by name and convert to bedpe format
echo -e "samtools view -bf 0x2 $outdir/${bam}_deduped_sorted.bam | \
bedtools bamtobed \
-i stdin \
-bedpe > $outdir/${bam}.bedpe\n" >> $shdir/${bam}.sh

### Format bedpe to bed
echo -e "Rscript $frags/motif_format_bedpe.R \
--id $bam \
--bedpe $outdir/${bam}.bedpe \
--outdir $outdir\n" >> $shdir/${bam}.sh

### Get FASTA sequences 
echo -e "bedtools getfasta \
-bedOut \
-fi $ref \
-bed $outdir/${bam}_5.bed > $outdir/${bam}_fasta_5.bed
bedtools getfasta \
-bedOut \
-fi $ref \
-bed $outdir/${bam}_3.bed > $outdir/${bam}_fasta_3.bed\n" >> $shdir/${bam}.sh

### Convert FASTA to motif context frequencies
echo -e "Rscript $frags/motif_get_contexts.R \
--id $bam \
--fasta_5 $outdir/${bam}_fasta_5.bed \
--fasta_3 $outdir/${bam}_fasta_3.bed \
--outdir $outdir\n" >> $shdir/${bam}.sh

### Remove intermediate files
echo -e "if [[ -f "$outdir/${bam}_raw.txt" ]]
then
echo -e \"Output completed sucessfully\"
rm $outdir/${bam}*.bam*
rm $outdir/${bam}_metrics.txt
rm $outdir/${bam}*.bed*
else
echo -e \"Errors in run\"
fi" >> $shdir/${bam}.sh
done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 10:00:00 -p all --mem 24G $file
done
