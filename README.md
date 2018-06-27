# Genomics_Pipeline
Genomic Pipelines (Mapping + Variant Calling) from NGS Data

## Getting Started
### Prerequisites
* UNIX/Linux
* Python 3 or higher 
* <a href="https://broadinstitute.github.io/picard">Picard</a>
* <a href="https://software.broadinstitute.org/gatk/documentation/quickstart">GenomeAnalysisToolkit (GATK 3.5-0)</a>
* <a href="https://sourceforge.net/projects/varscan/files/">VarScan.v2.3.9</a>

### Getting your clone
```
$git clone https://github.com/MBaysanLab/GenomicPipeline
```
If you do not have git you can download zipped Genomic Pipeline <a href="https://github.com/MBaysanLab/GenomicPipeline/archive/master.zip">here.</a>

### Usage

```
$cd genomics_pipeline
$python gui.py
```
You can now
* Select mapping, variant calling or full pipeline,
* Enter your working directory and decide which algorithm you want to use for mapping and variant calling,
* Choose a library ID and germline folder.
Then you can start processing your data via submit button. The final VCF file will be saved in your working directory.

### Contact