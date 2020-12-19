## Comparative analysis of alignment algorithms in mouse CD4-positive T cells 

### 1. Raw data 

- **GEO number**: [GSE128615](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128615)

- **SRA number**: [SRP188970](https://www.ncbi.nlm.nih.gov/sra?term=SRP188970)

- **Study design**: To determine the influence of differential Kdm6a expression in immune cells, whole transcriptome analysis for CD4+ T cells from WT and Kdm6a cKO mice were performed using RNA-Seq.

- **Reference**: [J Clin Invest. 2019 Aug 12;129(9):3852-3863. doi: 10.1172/JCI126250.](https://pubmed.ncbi.nlm.nih.gov/31403472)


#### 1-1. Raw data download 

- rawdata_download.sh

```bash
#!/bin/bash

# Numbers of SRA data:

# SRR8757083: WT-rep1
# SRR8757085: WT-rep2
# SRR8757086: WT-rep3
# SRR8757084: cKO-rep1
# SRR8757088: cKO-rep2
# SRR8757089: cKO-rep3


mkdir rawdata

cd rawdata

SRR_number=(3 4 5 6 8 9)

for x in ${SRR_number[*]}
do
    prefetch SRR875708$x

done


cd ..
```

#### 1-2. SRA to FASTQ files 

- sra_fastq.sh

```bash
#!/bin/bash


# Numbers of SRA data:

# SRR8757083: WT-rep1
# SRR8757085: WT-rep2
# SRR8757086: WT-rep3
# SRR8757084: cKO-rep1
# SRR8757088: cKO-rep2
# SRR8757089: cKO-rep3

cd rawdata

SRR_number=(3 4 5 6 8 9)


for x in ${SRR_number[*]}
do

    fastq-dump --split-files SRR875708$x/SRR875708$x.sra 


done


cd ..
```

#### 1-3. Renaming 

- name_change.sh

```bash
#!/bin/bash

# SRR8757083: WT-rep1
# SRR8757085: WT-rep2
# SRR8757086: WT-rep3
# SRR8757084: cKO-rep1
# SRR8757088: cKO-rep2
# SRR8757089: cKO-rep3

cd rawdata

mv SRR8757083_1.fastq WT-rep1.fastq
mv SRR8757085_1.fastq WT-rep2.fastq
mv SRR8757086_1.fastq WT-rep3.fastq
mv SRR8757084_1.fastq cKO-rep1.fastq
mv SRR8757088_1.fastq cKO-rep2.fastq
mv SRR8757089_1.fastq cKO-rep3.fastq

cd ..
```

#### 1-4. Gunzip

- fastq_gunzip.sh

```bash
#!/bin/bash



cd rawdata

gzip -v -k *.fastq 



cd ..
```

### 2. Conda environment 

- Reused from [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)

#### 2-1. Tools for alignment/mapping

- **HISAT2**: http://daehwankimlab.github.io/hisat2
- **Samtools**: http://www.htslib.org/doc/#manual-pages
- **STAR**: https://github.com/alexdobin/STAR
- **Salmon**: https://salmon.readthedocs.io/en/latest
- **bedtools**: https://bedtools.readthedocs.io/en/latest
- **gawk**: https://www.gnu.org/software/gawk/manual/gawk.html

#### 2-2. Tools for counting and differential expression (DE) analysis

- **DESeq2**: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- **Tximport**: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
- **Rsubread**: https://bioconductor.org/packages/release/bioc/html/Rsubread.html

### 3. Reference files 

- Same references as [this project](https://github.com/Mira0507/mouse_index_m25/blob/master/README.md)
- [GENCODE](https://www.gencodegenes.org/mouse) GRCm38.p6 m25

### 4. STAR alignment

#### 4-1. Indexing

- Reused from [this project](https://github.com/Mira0507/mouse_index_m25/blob/master/README.md)

#### 4-2. Alignment

- star_alignment.sh

```bash

#!/bin/bash

# Define directory names 
outdir=star_output   # directory storing output files
indir=../rawdata        # input directory (fastq files)
refdir=/home/mira/Documents/programming/Bioinformatics/mouse_reference     # reference directory (absolute path needed!)
indexdir=star_index    # index directory
genome=*genome.fa      # reference file
gtf=*.gtf              # GTF file

samples=(WT-rep{1..3} cKO-rep{1..3})

mkdir $outdir 

cd $outdir

for read in ${samples[*]} 
do 
    STAR --runThreadN 16 --runMode alignReads --genomeDir $refdir/$indexdir --readFilesCommand zcat --sjdbGTFfile $refdir/$gtf -sjdbOverhang 100 --readFilesIn $indir/${read}.fastq.gz --outFileNamePrefix $read --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic --chimOutType Junctions
done

cd ..
# --readFilesIn: For paired-end reads, use comma separated list for read1, followed by space, followed by comma separated list for read2
# --readFilesCommand: required when the input files are .gzip format (e.g. --readFilesCommand zcat, --readFilesCommand gunzip -c, or --readFilesCommand bunzip2)
# Note: the output files are generated in the current directory
```

### 5. Salmon

#### 5-1. Indexing 

- Reused from [this project](https://github.com/Mira0507/mouse_index_m25/blob/master/README.md)

#### 5-2. Mapping

- salmon_map.sh

```bash
#!/bin/bash



# Define index file directory
ind=~/Documents/programming/Bioinformatics/mouse_reference/salmon_index/gencode_index

out=salmon_output

# Define file names 
samples=(WT-rep{1..3} cKO-rep{1..3})

mkdir $out
cd $out

for read in ${samples[*]}

do
    salmon quant -i $ind -l A --gcBias --seqBias -r ~/Documents/programming/Bioinformatics/mouse_cd4/rawdata/${read}.fastq.gz -p 16 --validateMappings -o ${read}.salmon_quant
done

cd ..
```

### 6. HISAT2

#### 6-1. Indexing

- Reused from [this project](https://github.com/Mira0507/mouse_index_m25/blob/master/README.md)

#### 6-2. Alignment

- hisat2_align.sh

```bash

#!/bin/bash 

# Define directory and sample names
refdir=../mouse_reference/hisat2_index/index   # reference directory
samples=(WT-rep{1..3} cKO-rep{1..3})                     # sample names
outdir=hisat2_output                          # output directory
indir=rawdata                                 # input directory

mkdir $outdir 


for read in ${samples[*]} 

do 

    hisat2 -q -p 16 --seed 23 -x $refdir -U $indir/$read.fastq -S $outdir/$read.sam 

done 
```

#### 6-3. Converting SAM to BAM 

- Uses **Samtools**
- samtobam.sh

```bash
#!/bin/bash

cd hisat2_output

input=(WT-rep{1..3} cKO-rep{1..3})

for f in ${input[*]}

do
    samtools view -bS -@ 16 $f.sam > $f.bam 
done 


cd ..
```

#### 6-4. Sorting 

- Uses **Samtools**
- hisat2_sort.sh

```bash
#!/bin/bash

cd hisat2_output

input=(WT-rep{1..3} cKO-rep{1..3})  

for f in ${input[*]}

do
    samtools sort -@ 16 $f.bam -o $f.sorted.bam

done 


cd ..

# Delete SAM files after samtools run
```
