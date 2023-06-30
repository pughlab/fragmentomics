INPUTDIR=/cluster/projects/pughlab/hereditary/external_data/TGL49_CHARM/LFS/LFS_WG/bams
basedir=/cluster/projects/pughlab/hereditary/projects/CHARM/LFS/repeat_masker
shdir=$basedir/sh_scripts
outdir=$basedir/output
sitedir=/cluster/projects/pughlab/hereditary/projects/CHARM/repeat_masker
picard_dir=/cluster/tools/software/picard/2.10.9

mkdir -p $shdir
mkdir -p $outdir

cd $sitedir
ls *.bed > $shdir/sites

cd $INPUTDIR
ls *.bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
mv bam bams
sed 's/....$//' sites > site
mv site sites

for bam in $(cat bams); do

name=$(echo $bam | sed 's/\_WG.*/_WG/')
echo $bam
echo $name

mkdir -p $outdir/$name

echo -e "#!/bin/bash\n
module load picard\n" > $shdir/${name}.sh

echo -e "java -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
I=$INPUTDIR/${bam}.bam \
O=$outdir/$name/${name}_rm_genome.txt \
H=$outdir/$name/${name}_rm_genome.pdf \
M=0 \
W=600" >> $shdir/${name}.sh

for site in $(cat sites); do
echo -e "### $site ###
#samtools view -bh -L $sitedir/${site}.bed $INPUTDIR/${bam}.bam > $outdir/$name/${site}.bam
#samtools index $outdir/$name/${site}.bam" >> $shdir/${name}.sh

echo -e "#java -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
I=$outdir/$name/${site}.bam \
O=$outdir/$name/${name}_rm_${site}.txt \
H=$outdir/$name/${name}_rm_${site}.pdf \
M=0 \
W=600" >> $shdir/${name}.sh

echo -e "#rm $outdir/$name/${site}.bam*\n" >> $shdir/${name}.sh

done
done

cd $shdir

ls *.sh > files
for file in $(cat files); do
sbatch -c 1 -t 24:00:00 -p all --mem 16G $file
done
