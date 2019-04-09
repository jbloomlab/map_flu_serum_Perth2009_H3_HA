# Uploading FASTQ files to the SRA
The script instructions and for how the raw sequencing files (FASTQ files) were uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra):

1. Go to the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) and manually create a *BioProject* and appropriate *BioSamples*.
   These are then entered into the Python script [make_sra_forms.py](make_sra_forms.py) near the top to specify the *BioProject* and *BioSamples*.

2. Run [make_sra_forms.py](make_sra_forms.py), a Python script that reads in the [../data/serum_info.yaml](../data/serum_info.yaml) and [../data/sample_list.csv](../data/sample_list.csv) files and uses them to generate the files needed for submission.
   Specifically, it generates the submission form [submissionform_2019_Perth_serum_MAP.tsv](submissionform_2019_Perth_serum_MAP.tsv) and the FASTQ `tar` file [2019_Perth_serum_MAP.tar](2019_Perth_serum_MAP.tar).
   The `tar` file takes a long time to generate, so you might first want to run the script with the `make_tar` option set to `False` to make sure it is working.

3. Once you have generated the submission form and `tar` of the FASTQ files, **Juhye, what happens next?**
   


