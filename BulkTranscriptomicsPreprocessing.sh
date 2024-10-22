# Bulk transcriptomics preprocessing 

# FastQC
~/Softwares/FastQC/fastqc -o ~/FASTQC -f fastq --no-extract -t 64 ~/FASTQ/*.fq.gz

# Genome Generation using STAR 
cd ~/Genome
~/Softwares/STAR/STAR-2.7.11a/bin/MacOSX_x86_64/STAR --runThreadN 64 --runMode genomeGenerate --genomeDir ~/Genome --genomeFastaFiles ~/Genome/GRCm39.genome.fa --sjdbGTFfile ~/Genome/gencode.vM33.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 149

# STAR Alignment
for infile in *_1.fq.gz
do 
base=$(basename ${infile} _1.fq.gz)
~/Softwares/STAR/STAR-2.7.11a/bin/MacOSX_x86_64/STAR --runThreadN 64 --genomeDir ~/Genome --twopassMode Basic --outSAMstrandField intronMotif --sjdbOverhang 149 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn ~/FASTQ/${infile} ${base}_2.fq.gz --outFileNamePrefix ~/STAR/${base} --outTmpDir ~/tmp
done

# GTF file conversion to BED file for Infer-Experiment
cd ~/BAM
/usr/local/bin/gtf2bed < ~/Genome/gencode.vM33.chr_patch_hapl_scaff.annotation.gtf > Mus32.bed

# InferExperiment to check for strandedness
cd ~/BAM
infer_experiment.py -r Mus32.bed -i C57_1FAligned.sortedByCoord.out.bam

# Index bam file to visualise the mapping results using IGV
cd ~/BAM
for infile in *.bam
do
~/Softwares/samtools/samtools-1.16.1/samtools index ${infile}
done

# FeatureCounts
cd ~/BAM
~/Softwares/subread-2.0.3-macOS-x86_64/bin/featureCounts -p --countReadPairs -t exon -g gene_id -s 2 -M -a ~/Counts/gencode.vM34.chr_patch_hapl_scaff.annotation.gtf -o ~/Counts/RawCounts_Full.txt #read in all bam files following this code by listing

# StringTie
cd ~/BAM
for infile in *.bam
do
base=$(basename ${infile} .bam)
~/Softwares/stringtie/stringtie -A ~/StringTie/${base}.txt -G ~/StringTie/gencode.vM34.chr_patch_hapl_scaff.annotation.gtf ${infile}
done 


