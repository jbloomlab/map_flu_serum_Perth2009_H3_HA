# Mapping of anti-flu serum against the Perth/2009 H3 HA
Mutational antigenic profiling of Perth/2009 H3 HA codon mutant libraries against ferret and human sera.

Study by Juhye Lee, Rachel Eguia, and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html).

## Quick summary
Look at the [notebook results](results/analysis_notebook.md) for an overview of the results.

## Running the analysis
The main analysis is performed by the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb).
The Markdown results of the notebook are shown in [results/analysis_notebook.md](results/analysis_notebook.md).

To run [analysis_notebook.ipynb](analysis_notebook.ipynb) and generate the [Markdown results](results/analysis_notebook.md), run the bash script [run_notebook.bash](run_notebook.bash) with:

    ./run_notebook.bash
    
On the Hutch cluster, you probably want to submit this job using [slurm](https://slurm.schedmd.com/), which you can simply do with:

    sbatch -t 2-0 -p largenode -c 16 --mem=100000 run_notebook.bash

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
      - *R1*: glob pattern to R1 FASTQ files on Hutch server
      - **to be added**: eventually we will also add the SRA accesssion when the FASTQ files are uploaded to the SRA so there can be an option to run the analysis on either the FASTQ files on Hutch server or those downloaded from SRA.
  
  - [data/Perth09_HA_reference.fa](data/Perth09_HA_reference.fa):
    FASTA file giving the sequence of the wildtype Perth/2009 HA used in the experiments.
  
  - [data/H3renumbering_scheme.csv](data/H3renumbering_scheme.csv):
    A CSV file that maps sequential 1, 2, ... numbering of the Perth/2009 HA protein sequence (*original* column) to the standard H3 HA numbering scheme (*new* column).
  
## Results
Results are placed in the [./results/](results) subdirectory.
Many of the results files are **not** tracked in this GitHub repo since they are very large.
However, the following results are tracked:

  - [results/analysis_notebook.md](results/analysis_notebook.md): Markdown rendering of results of the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb)

## Other subdirectories
Other subdirectories have information relevant to the study that are not part of the main pipeline described above:

 - [./neutralization_assays/](neutralization_assays) has data and analysis relevant to the neutralization assays.
   See [./neutralization_assays/README.md](neutralization_assays/README.md) for details.

 - [./SRA_upload/](SRA_upload) has information on how the sequencing data are uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).
   See [./SRA_upload/README.md](SRA_upload/README.md) for details.
