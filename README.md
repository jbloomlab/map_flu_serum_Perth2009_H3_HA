# Mapping of anti-flu serum against the Perth/2009 H3 HA
Mutational antigenic profiling of Perth/2009 H3 HA codon mutant libraries against ferret and human sera.

Study by Juhye Lee, Rachel Eguia, and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html).

## Quick summary
- Analysis of mutational antigenic profiling data: [results/notebooks/analyze_map.md](results/notebooks/analyze_map.md)

- Analysis of neutralization curves: [results/notebooks/analyze_neut.md](results/notebooks/analyze_neut.md)

## Running the analysis
The main analysis is performed primarily by a series of [Jupyter notebooks](https://jupyter.org/) and Python scripts:

  1. [analyze_map.ipynb](analyze_map.ipynb): analyzes mutational antigenic profiling

  2. [analyze_neut.ipynb](analyze_neut.ipynb): analyzes neutralization assays

  3. [parameterize_map_on_struct.py](parameterize_map_on_struct.py): parameterizes the template Jupyter notebook [map_on_struct_template.ipynb](map_on_struct_template.ipynb) to show structures for each type of sera.

To run the three steps above, execute the bash script [run.bash](run.bash) with:

    ./run.bash
    
On the Hutch cluster, you can also submit this script using [slurm](https://slurm.schedmd.com/) with:

    sbatch -p largenode -c 16 --mem=100000 run.bash

Running these three steps creates Jupyter notebooks that map the immune selection onto the structure using [dms_struct](https://jbloomlab.github.io/dms_struct) (which is a wrapper around [nglview](https://github.com/arose/nglview)).
These notebooks are in [results/notebooks/structs](results/notebooks/structs).
They can be used to render interactive views of the structures via [mybinder](https://mybinder.org/) with [appmode](https://github.com/oschuett/appmode) by clicking on the links in the *Quick summary* section above.
You should also run them interactively cell-by-cell (giving time for each structure to render) in order to create static structure images.

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

  - [data/H3_site_to_PDB_4o5n.csv](data/H3_site_to_PDB_4o5n.csv):
    A CSV file that matches the H3 HA numbering to the site numbers and PDB chains in PDB structure [4o5n](https://www.rcsb.org/structure/4O5N).

  - [data/neut_assays](data/neut_assays):
    Data from neutralization assays:
      - [data/neut_assays/neut_config.yaml](data/neut_assays/neut_config.yaml): Details on neutralization assays for each serum, with Excel path relative to top-level analysis directory.
      - [data/plate_reader_data/](data/plate_reader_data/): The plate reader data (Excel format).

  - [data/figure_config.yaml](data/figure_config.yaml):
    YAML file giving specifications for fine-tuned figures showing logo plot zooms and neutralization curves. 
  
## Results
Results are placed in the [./results/](results) subdirectory.
Many of the results files are **not** tracked in this GitHub repo since they are very large.
However, the following results are tracked:

  - [results/notebooks/analyze_map.md](results/notebooks/analyze_map.md): Markdown rendering of the notebook analyzing the mutational antigenic profiling

  - [results/notebooks/analyze_neut.md](results/notebooks/analyze_neut.md): Markdown rendering of the notebook analyzing the neutralization assays

  - [results/notebooks/structs](results/notebooks/structs): Jupyter notebooks that render the immune selection onto interactive structure widgets.

  - [results/renumbered_codoncounts](results/renumbered_codoncounts): the counts of each mutation in each sample, after renumbering to H3 numbering.

## Other subdirectories
Other subdirectories in the repo are:

 - [./SRA_upload/](SRA_upload) has information on how the sequencing data were uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).
   See [./SRA_upload/README.md](SRA_upload/README.md) for details.

 - [./binder](binder) has the configuration for visualizing the protein structure notebooks on [mybinder](https://mybinder.org/).
