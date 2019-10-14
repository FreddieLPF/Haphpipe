# Haphpipe User Guide

---
---

## Basics:

*Made by:* Uzma Rentia

*My virus is:* HCV

*My links are:*

*Associated paper pararaphs.*
*Paragraph about study:*
Chronic hepatitis C virus (HCV) infection is the most common cause of end-stage liver disease leading to liver transplantation in the country, with recurrent HCV disease being more aggressive during the course of immunosuppression following transplants. MBL-HCV1, a monoclonal antibody that binds to the highly-conserved epitope of the HCV E2 envelope glycoprotein, has been show to suppress viral load in chimpanzee models. A viral load uptick following treatment is concurrent with the emergence of antibody-resistance virus and is most consistent with amino acid alterations at positions and 417 (N417S) within the MBL-HCV1 epitope. This study aimed to track the development of post-transplant HCV variants with and without MBL-HCV1. As in the chimp studies, it found that MBL-HCV1 was successful until the development of the aforementioned variants, which was confirmed by Sanger sequencing. Furthermore, high-throughput sequencing sensitive enough to detect MAb-resistance associated variants before transplantation confirmed that the resistant variants of note were not present prior to treatment. 

*Paragraph abut NGS technique and data:*

## Questions to answer in your own words
Do not feel like you need a "perfect answer" - that is not what we're going for. Explain things in your own words and as if you're explaining them to your peers or another undergrad you "tutor". We're trying to make this approachable not perfect, scientific explainations. Also, you don't have to answer these *prior* to starting, just eventually or throughout your work on this.

1. What is a directory structure? Please explain and create a diagram / picture to accompany your explaination.
2. What is the difference between next-generation sequencing and Sanger sequencing? The two are very similar except that NGS can sequence multiple DNA fragments at once. Additionally, NGS is more sensitive to lower volumes of DNA
3. What is the difference between the two pipelines, hp01 and hp02? Explain what each step of the pipelines accomplishes and why that step is necessary.
4. Did you need to learn how to bash script for this? What do you feel like you weren't prepared for after reading the introduction part we provided you with?



## Installing HAPHPIPE

Insert instructions in your own words and code.

## Activating HAPHPIPE

Insert how you activate HAPHPIPE and use it in your own words and code.


## HAPHPIPE suite
Please explain what the haphpipe suite is. What is it's purpose how do you use it? What is it good for?

For each of the stages below, please explain the stage. List the files needed for the command, what the command does. The options for the command. An example of how to execute the command.

### hp_reads
        sample_reads
        trim_reads
        join_reads
        ec_reads
### hp_assemble
        assemble_denovo
        assemble_amplicons
        assemble_scaffold
        align_reads
        call_variants
        vcf_to_consensus
        refine_assembly
        finalize_assembly
### hp_haplotype
        predict_haplo
        ph_parser
### hp_annotate
        pairwise_align
        extract_pairwise
        annotate_from_ref
        
## Example usage
Your virus with both pipelines. Document all code and explanation.  

## Helpful resources
List all helpful topics you think or did address here with links / explainations you felt were helpful. Use bullet points.

- introduction to NGS
	- https://learn.gencore.bio.nyu.edu/variant-calling/
- inroduction to bash
	
## FAQ
List all questions you had and provide links / explainations you felt were helpful to answer these questions. Use topics and numbering.
