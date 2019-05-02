# Mapping of anti-flu serum against the Perth/2009 H3 HA
Mutational antigenic profiling of Perth/2009 H3 HA codon mutant libraries against ferret and human sera.

Study by Juhye Lee, Rachel Eguia, and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html).

## Quick summary
Look at the [Markdown notebook results](results/notebooks) for an overview of the results:

  - Analysis of mutational antigenic profiling data: [results/notebooks/analyze_map.md](results/notebooks/analyze_map.md)

  - Analysis of neutralization: [results/notebooks/analyze_neut.md](results/notebooks/analyze_neut.md)

## Running the analysis
The main analysis is performed by the following [Jupyter notebooks](https://jupyter.org/):

  - [analyze_map.ipynb](analyze_map.ipynb): analyzes mutational antigenic profiling

  - [analyze_neut.ipynb](analyze_neut.ipynb): analyzes neutralization assays

The easiest way to look at the results is to view the Markdown output of the notebooks as described above.

To run the notebooks and generate the Markdown output, run the bash script [run_nbs.bash](run_nbs.bash) with:

    ./run_nbs.bash
    
On the Hutch cluster, you probably want to submit this job using [slurm](https://slurm.schedmd.com/), which you can simply do with:

    sbatch -p largenode -c 16 --mem=100000 run_nbs.bash

## Configuring the analysis
The configuration for the analysis is in a separate file, [config.yaml](config.yaml). 
This file defines key variables for the analysis, and should be self-explanatory. 
The [config.yaml](config.yaml) file points to several files in the [./data/](data) subdirectory that specify essential data for the analysis:

  - [data/serum_info.yaml](data/serum_info.yaml):
    YAML file that gives information on all of the serum samples used for selections.
    For each serum there is an entry with the label used in the experiments, then:
      - *name*: a more informative name used when displaying results
      - *description*: description of the serum
      - *group*: group of samples to which serum belong
      - *species*: species from which serum is derived (if relevant)
      - *vaccination*: information of vaccination status (if relevant)

  - [data/sample_list.csv](data/sample_list.csv):
    CSV file giving each sample that was deep sequenced.
    Columns are:
      - *sample*: sample label used in experiments
      - *serum*: serum used for selection
      - *library*: viral library, using simple 1, 2, 3 naming rather than the more confusing library codes used to label experiments
      - *date*: day when sequencing was done
      - *serum_dilution*: dilution of serum used; this includes the 1:4 dilution used during the RDE treatment of the serum. For antibodies, it is the concentration in ug/ml.
      - *percent_infectivity*: percent of viral library retaining infectivity
      - *R1*: glob pattern to R1 FASTQ files on Hutch server; the R2 file names are guessed from the R1 names. If [config.yaml](config.yaml) sets *seq_data_source* to *R1* then there must be a valid R1 file glob for all samples; otherwise this column is ignored.
      - *SRA_accession*: the accession number on the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA) for the sequencing for this sample. If [config.yaml](config.yaml) sets *seq_data_source* to *SRA_accession* then there must be a valid accession for all samples; otherwise this column is ignored.
  
  - [data/Perth09_HA_reference.fa](data/Perth09_HA_reference.fa):
    FASTA file giving the sequence of the wildtype Perth/2009 HA used in the experiments.
  
  - [data/H3renumbering_scheme.csv](data/H3renumbering_scheme.csv):
    A CSV file that maps sequential 1, 2, ... numbering of the Perth/2009 HA protein sequence (*original* column) to the standard H3 HA numbering scheme (*new* column).

  - [data/neut_assays](data/neut_assays):
    Data from neutralization assays.
    Specifically:
      - [data/neut_assays/neut_config.yaml](data/neut_assays/neut_config.yaml): Details on neutralization assays for each serum, with Excel path relative to top-level analysis directory.
      - [data/plate_reader_data/](data/plate_reader_data/): The plate reader data (Excel format).

  - [data/mutation_colors_and_markers.yaml](data/mutation_colors_and_markers.yaml):
    YAML file giving specifications for colored logo plot zooms and neutralization assay plots. 
  
## Results
Results are placed in the [./results/](results) subdirectory.
Many of the results files are **not** tracked in this GitHub repo since they are very large.
However, the following results are tracked:

  - [results/notebooks/analyze_map.md](results/notebooks/analyze_map.md): analysis of mutational antigenic profiling

  - [results/notebooks/analyze_neut.md](results/notebooks/analyze_neut.md): analysis of neutralization assays

## Other subdirectories
Other subdirectories have information relevant to the study that are not part of the main pipeline described above:

 - [./SRA_upload/](SRA_upload) has information on how the sequencing data were uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).
   See [./SRA_upload/README.md](SRA_upload/README.md) for details.
