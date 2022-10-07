INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
basedir=/cluster/projects/pughlab/projects/CHARM/LFS/nucleosome_peaks
ref=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
picard_dir=/cluster/tools/software/picard/2.10.9
frags=/cluster/projects/pughlab/bin/fragmentomics
peaks=/cluster/projects/pughlab/projects/CHARM/nucleosome_peaks/hg38
outdir=$basedir/output/fragments
shdir=$basedir/sh_scripts/fragments

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

### Subset the bam file for insert sizes of 167bp
echo -e "samtools view -h $INPUTDIR/${bam}.bam | awk 'substr(\$0,1,1)==\"@\" || (\$9==167) || (\$9==-167)' | \
samtools view -b > $outdir/${bam}.bam
samtools index $outdir/${bam}.bam\n" >> $shdir/${bam}.sh

### Remove duplicates
echo -e "java -jar $picard_dir/picard.jar MarkDuplicates \
I=$outdir/${bam}.bam \
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

### Split bedpe into chromosomes
echo -e "$frags/splitBedpe.sh $outdir/${bam}.bedpe\n" >> $shdir/${bam}.sh

### Calculate the distance from closest peak
echo -e "Rscript $frags/R/nucleosome_peaks_distance.R \
--id $bam \
--path $outdir \
--peaks $peaks \
--outdir $outdir\n" >> $shdir/${bam}.sh

### Remove intermediate files
echo -e "if [[ -f "$outdir/${bam}_peak_distance.txt" ]]
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
