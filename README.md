
# Haphpipe User Guide

### Basics:

*Made by: Freddie Pengfei Li*

*My virus is: HIV(human immunodeficiency viruses)*

*My links are: [MyGithub](https://github.com/FreddieLPF/haphpipe_bioinformaticsclass)*

*
Paragraph abut NGS technique and data: 
Advantages of NGS include: 
1. Higher sensitivity to detect low-frequency variants.
2. Faster turnaround time for high sample volumes.
3. Comprehensive genomic coverage Lower limit of detection.
4. Higher throughput with sample multiplexing.
5. Ability to sequence hundreds to thousands of genes or gene regions simultaneously
*

  

## Core Questions

**1. What is a directory structure? Please explain and create a diagram / picture to accompany your explanation.**

In computing, a directory structure is the way an operating system's file system and its files are displayed to the user.

As researchers and software engineers, we should make sure that someone who is unfamiliar with our project is able to look at our computer files and understand in detail what we did and why.

We all have the rough experience where after a few months, we may simply not remember what we were up to when we created a particular set of files, or we may be forgetful about what conclusions were drew. We will either have to then spend time reconstructing previous experiments or lose whatever insights gained from those experiments.

So maintaining a well-structured directory is essential to both software development and bioinformatics.


 **2. What is the difference between next-generation sequencing and Sanger sequencing?**
  
The critical difference between Sanger sequencing and NGS is sequencing volume. While the Sanger method only sequences a single DNA fragment at a time, NGS is massively parallel, sequencing millions of fragments simultaneously per run. This high-throughput process translates into sequencing hundreds to thousands of genes at one time. NGS also offers greater discovery power to detect novel or rare variants with deep sequencing.

  

**3. What is the difference between the two pipelines, hp01 and hp02? Explain what each step of the pipelines accomplishes and why that step is necessary.**

- hp_assemble_01

>This pipeline implements denovo assembly. Reads are first trimmed (trim_reads) and used as input for denovo assembly (assemble_denovo). The denovo assembly stage automatically performs error detection, but trimmed reads are also error corrected (ec_reads). The assembled contigs are used as input for amplicon assembly (assemble_amplicons) along with reference FASTA and GTF files. The assembly is then iteratively refined up to five times (refine_assembly) by mapping corrected reads to the assembled FASTA file and lastly finalized (finalize_assembly), resulting in a FASTA file with final consensus sequences, final VCF, and aligned BAM file.

  

- hp_assemble_02  
>This pipeline implements reference-based mapping assembly. Reads are first trimmed (trim_reads) and error-corrected (ec_reads). The corrected reads are used as input for reference-based mapping assembly (refine_assembly) for up to five iterations. Lastly, the assembly is finalized (finalize_assembly) by mapping reads onto the refined reference sequence. The final output is a FASTA file with final consensus sequences, final VCF, and aligned BAM file.

**4. Did you need to learn how to bash script for this? What do you feel like you weren't prepared for after reading the introduction part we provided you with?**

  Some basic bash-scripting knowledge is required to run haphpipe. Personally as a computer science major, I think it would be better if I was given more background about the biology part of the pipeline, namely: What does the terminology mean? How does the pipeline run? How does the script apply to research and science field? What result would we get after running?

  
## Installing HAPHPIPE

**1.  Conda**
Conda is an open source package management system and environment management system that runs on Windows, macOS and Linux. Conda as a package manager helps you find and install packages. Do conda info -h whenever you are not sure and a manual will pop up.

**2. Conda channels**
Channels are the locations where packages are stored. They serve as the base for hosting and managing packages. Conda packages are downloaded from remote channels, which are URLs to directories containing Conda packages. So we need to add channels first: 
```
conda config --add channels channel_name
```
_channel_name_ would be 

```
- defaults 
- bioconda
- conda-forge
 ```

**3. Create a conda environment with the necessary dependencies:**
```
conda create -n haphpipe_freddie \
 python \
 future \
 pyyaml \
 biopython \
 seqtk \
 bowtie2 \
 bwa \
 flash \
 freebayes \
 mummer \
 picard \
 trimmomatic \
 samtools=1.9 \
 gatk=3.8 \
 spades \
 blast \
 sierrapy
 ```

**_Note_**:  haphpipe_freddie is an example, you can rename the name of your Conda environment.

when ```Proceed ([y]/n)?``` shows up, type ```y``` when you are sure you want to proceed and install.

  
**5.Install GATK:**

Due to license restrictions, bioconda cannot distribute and install GATK directly. To fully install GATK, you need to download a licensed copy of GATK (version 3.8-0) from the [Broad Institute](https://software.broadinstitute.org/gatk/download/archive).

And then register the package using gatk3-register:

```gatk3-register /path/to/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2```

This will copy GATK into your conda environment.

Note: HAPHPIPE was developed and tested using GATK 3.8.

  **6. Install Haphpipe:**

```
pip install git+git://github.com/gwcbi/haphpipe.git
```

Now you are all set, ```Successfully installed haphpipe``` will pop up on your screen.

## Activating HAPHPIPE

To activate this environment, use

```conda activate haphpipe_freddie```

To deactivate an active environment, use

``` conda deactivate```

  

## HAPHPIPE suite

### Stages
Each stage can be run on its own. Stages are grouped into 4 categories: hp_reads, hp_assemble, hp_haplotype, and hp_annotate.

More detailed description of command line options for each stage are available in the [wiki](https://github.com/gwcbi/haphpipe/wiki). 

To view all available stages in HAPHPIPE, run:

```
haphpipe -h

```

  

  

  

### Reads

  

Stages to manipulate reads and perform quality control. Input is reads in FASTQ format, output is modified reads in FASTQ format.

  

##### sample_reads

  

Subsample reads using seqtk ([documentation](https://github.com/lh3/seqtk)). Input is reads in FASTQ format. Output is sampled reads in FASTQ format.

Example to execute:

```

haphpipe sample_reads --fq1 read_1.fastq --fq2 read_2.fastq --nreads 1000 --seed 1234

```

  

##### trim_reads

  

Trim reads using Trimmomatic ([documentation](http://www.usadellab.org/cms/?page=trimmomatic)). Input is reads in FASTQ format. Output is trimmed reads in FASTQ format.

Example to execute:

```

haphpipe trim_reads --fq1 read_1.fastq --fq2 read_2.fastq

```

  

##### join_reads

  

Join reads using FLASH ([paper](https://www.ncbi.nlm.nih.gov/pubmed/21903629)). Input is reads in FASTQ format. Output is joined reads in FASTQ format.

Example to execute:

```

haphpipe join_reads --fq1 trimmed_1.fastq --fq2 trimmed_2.fastq

```

  

##### ec_reads

  

Error correction using SPAdes ([documentation](http://cab.spbu.ru/software/spades/)). Input is reads in FASTQ format. Output is error-corrected reads in FASTQ format.

Example to execute:

```

haphpipe ec_reads --fq1 trimmed_1.fastq --fq2 trimmed_2.fastq

```

  

### Assemble

  

Assemble consensus sequence(s). Input reads (in FASTQ format) are assembled

using either denovo assembly or reference-based alignment.

Resulting consensus can be further refined.

  

##### assemble_denovo

  

Assemble reads via de novo assembly using SPAdes ([documentation](http://cab.spbu.ru/software/spades/)). Input is reads in FASTQ format. Output is contigs in FNA format.

Example to execute:

```

haphpipe assemble_denovo --fq1 corrected_1.fastq --fq2 corrected_2.fastq --outdir denovo_assembly --no_error_correction TRUE

```

##### assemble_amplicons

  

Assemble contigs from de novo assembly using both a reference sequence and amplicon regions with MUMMER 3+ ([documentation](http://mummer.sourceforge.net/manual/)). Input is contigs and reference sequence in FASTA format and amplicon regions in GTF format.

Example to execute:

```

haphpipe assemble_amplicons --contigs_fa denovo_contigs.fa --ref_fa refSequence.fasta --ref_gtf refAmplicons.gtf

```

  

##### assemble_scaffold

  

Scaffold contigs against a reference sequence with MUMMER 3+ ([documentation](http://mummer.sourceforge.net/manual/)). Input is contigs in FASTA format and reference sequence in FASTA format. Output is scaffold assembly, alligned scaffold, imputed scaffold, and padded scaffold in FASTA format.

Example to execute:

```

haphpipe assemble_scaffold --contigs_fa denovo_contigs.fa --ref_fa refSequence.fasta

```

  

##### align_reads

  

Map reads to reference sequence (instead of running de novo assembly) using Bowtie2 ([documentation](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)) and Picard ([documentation](https://broadinstitute.github.io/picard/)). Input is reads in FASTQ format and reference sequence in FASTA format.

Example to execute:

```

haphpipe align_reads --fq1 corrected_1.fastq --fq2 corrected _2.fastq --ref_fa refSequence.fasta

```

  

##### call_variants

  

Variant calling from alignment using GATK ([documentation](https://software.broadinstitute.org/gatk/download/archive)). Input is alignment file in BAM format and reference sequence in FASTA format (either reference from reference-based assembly or consensus final sequence from de novo assembly). Output is a Variant Call File (VCF) format file.

Example to execute:

```

haphpipe call_variants --aln_bam alignment.bam --ref_fa refSequence.fasta

```

  

##### vcf_to_consensus

  

Generate a consensus sequence from a VCF file. Input is a VCF file. Output is the consensus sequence in FASTA format.

Example to execute:

```

haphpipe vcf_to_consensus --vcf variants.vcf

```

  

##### refine_assembly

  

Map reads to a denovo assembly or reference alignment. Assembly or alignment is iteratively updated. Input is reads in FASTQ format and reference sequence (assembly or reference alignment) in FASTA format. Output is refined assembly in FASTA format.

Example to execute:

```

haphpipe refine_assembly --fq_1 corrected_1.fastq --fq2 corrected_2.fastq --ref_fa refSequence.fasta

```

  

##### finalize_assembly

  

Finalize consensus, map reads to consensus, and call variants. Input is reads in FASTQ format and reference sequence in FASTA format. Output is finalized reference sequence, alignment, and variants (in FASTA, BAM, and VCF formats, respectively).

```

haphpipe finalize_assembly --fq_1 corrected_1.fastq --fq2 corrected_2.fastq --ref_fa refined.fna

```

  

### Haplotype

  

Haplotype assembly stages. HAPHPIPE implements PredictHaplo ([paper](https://www.ncbi.nlm.nih.gov/pubmed/26355517)), although other haplotype reconstruction programs can be utilized outside of HAPHPIPE using the final output of HAPHPIPE, typically with the final consensus sequence (FASTA) file, reads (raw, trimmed, and/or corrected), and/or final alignment (BAM) file as input.

  

##### predict_haplo

  

Haplotype identification with PredictHaplo. Input is reads in FASTQ format and and reference sequence in FASTA format. Output is the longest global haplotype file and corresponding HTML file. __Note: PredictHaplo must be installed separately before running this stage.__

Example to execute:

```

haphpipe predict_haplo corrected_1.fastq --fq2 corrected_2.fastq --ref_fa final.fna

```

  

##### ph_parser

  

Return PredictHaplo output as a correctly formatted FASTA file. Input is the output file from __predict_haplo__ (longest global .fas file). Output is a correctly formatted FASTA file.

Example to execute:

```

haphpipe ph_parser best.fas

```

  

### Annotate

  

Stages to annotate consensus sequences.

  

##### pairwise_align

  

Apply correct coordinate system to final sequence(s) to facilitate downstream analyses. Input is the final sequence file in FASTA format, a reference sequence in FASTA format, and a reference GFT file. Output is a JSON file to be used in __extract_pairwise__.

Example to execute:

```

haphpipe pairwise_align --amplicons_fa final.fna --ref_fa refSequence.fasta --ref_gtf referenceSeq.gtf

```

  

##### extract_pairwise

  

Extract sequence regions from the pairwise alignment produced in __pairwise_align__. Input is the JSON file from __pairwise_align__. Output is either an unaligned nucleotide FASTA file, an aligned nucleotide FASTA file, an amino acid FASTA file, an amplicon GTF file, or a tab-separated values (TSV) file (default: nucleotide FASTA with regions of interest from GTF file used in __pairwise_align__).

Example to execute:

```

haphpipe extract_pairwise --align_json pairwise_aligned.json --refreg HIV_B.K03455.HXB2:2085-5096

```

  

##### annotate_from_ref

  

Annotate consensus sequence from reference annotation. Input is JSON file from __pairwise_align__ and reference GTF file.

Example to execute:

```

haphpipe annotate_from_ref airwise_aligned.json --ref_gtf referenceSeq.gtf

```

  

## Example usage

```
haphpipe_assemble_01 read1.fq.gz read2.fq.gz ../refs/HIV_B.K03455.HXB2.fasta ../refs/HIV_B.K03455.HXB2.gtf sampleID
```

```
haphpipe_assemble_02 read1.fq.gz read2.fq.gz ../refs/HIV_B.K03455.HXB2.amplicons.fasta sampleID
```

  

## Helpful resources

List all helpful topics you think or did address here with links / explanations you felt were helpful. Use bullet points.


- [introduction to NGS](https://learn.gencore.bio.nyu.edu/variant-calling)
- [introduction to bash](https://linuxconfig.org/bash-scripting-tutorial-for-beginners)
- [intro to github](https://github.com/gwcbi/HPC/blob/master/github.md)
- [markdown tutorial](https://www.markdowntutorial.com)
- [tips for command line](https://github.com/gwcbi/HPC/blob/master/commandline.md)
- [how to run interactive jobs (do this for all haphpipe modules and pipelines!](https://github.com/gwcbi/HPC/blob/master/interactive_jobs.md)
  

## FAQ

- [what is conda?](https://conda.io/en/latest/)

- [A review of bioinformatic pipeline frameworks](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429012/)

- [A Quick Guide to Organizing Computational Biology Projects](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

- [Directory structure](https://en.wikipedia.org/wiki/Directory_structure)

- [ngs-vs-sanger-sequencing](https://www.illumina.com/science/technology/next-generation-sequencing/ngs-vs-sanger-sequencing.html)
