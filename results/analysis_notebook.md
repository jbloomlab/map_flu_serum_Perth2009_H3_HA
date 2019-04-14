
<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Analysis-of-serum-mutational-antigenic-profiling" data-toc-modified-id="Analysis-of-serum-mutational-antigenic-profiling-1">Analysis of serum mutational antigenic profiling</a></span><ul class="toc-item"><li><span><a href="#Configure-analysis" data-toc-modified-id="Configure-analysis-1.1">Configure analysis</a></span><ul class="toc-item"><li><span><a href="#Import-Python-modules-/-packages" data-toc-modified-id="Import-Python-modules-/-packages-1.1.1">Import Python modules / packages</a></span></li><li><span><a href="#Get-config-variables" data-toc-modified-id="Get-config-variables-1.1.2">Get config variables</a></span></li><li><span><a href="#Get-sera-information" data-toc-modified-id="Get-sera-information-1.1.3">Get sera information</a></span></li><li><span><a href="#Get-list-of-samples" data-toc-modified-id="Get-list-of-samples-1.1.4">Get list of samples</a></span></li></ul></li><li><span><a href="#Download-deep-sequencing-data-if-needed" data-toc-modified-id="Download-deep-sequencing-data-if-needed-1.2">Download deep sequencing data if needed</a></span></li><li><span><a href="#Align-sequencing-and-count-mutations" data-toc-modified-id="Align-sequencing-and-count-mutations-1.3">Align sequencing and count mutations</a></span><ul class="toc-item"><li><span><a href="#Run-dms2_batch_bcsubamp" data-toc-modified-id="Run-dms2_batch_bcsubamp-1.3.1">Run <code>dms2_batch_bcsubamp</code></a></span></li><li><span><a href="#Plot-sequencing-and-mutation-counts-summaries" data-toc-modified-id="Plot-sequencing-and-mutation-counts-summaries-1.3.2">Plot sequencing and mutation counts summaries</a></span><ul class="toc-item"><li><span><a href="#Reads-and-barcodes-per-sample" data-toc-modified-id="Reads-and-barcodes-per-sample-1.3.2.1">Reads and barcodes per sample</a></span></li><li><span><a href="#Coverage-across-gene" data-toc-modified-id="Coverage-across-gene-1.3.2.2">Coverage across gene</a></span></li><li><span><a href="#Mutation-frequencies" data-toc-modified-id="Mutation-frequencies-1.3.2.3">Mutation frequencies</a></span></li><li><span><a href="#Check-for-oxidative-damage" data-toc-modified-id="Check-for-oxidative-damage-1.3.2.4">Check for oxidative damage</a></span></li></ul></li><li><span><a href="#Renumber-to-standard-HA-numbering" data-toc-modified-id="Renumber-to-standard-HA-numbering-1.3.3">Renumber to standard HA numbering</a></span></li></ul></li><li><span><a href="#Compute-immune-selection-on-mutations" data-toc-modified-id="Compute-immune-selection-on-mutations-1.4">Compute immune selection on mutations</a></span><ul class="toc-item"><li><span><a href="#Samples-to-compare-for-each-selection" data-toc-modified-id="Samples-to-compare-for-each-selection-1.4.1">Samples to compare for each selection</a></span></li><li><span><a href="#Choose-measure-of-immune-selection" data-toc-modified-id="Choose-measure-of-immune-selection-1.4.2">Choose measure of immune selection</a></span></li><li><span><a href="#Compute-immune-selection" data-toc-modified-id="Compute-immune-selection-1.4.3">Compute immune selection</a></span></li><li><span><a href="#Get-all-selection-information-in-one-data-frame" data-toc-modified-id="Get-all-selection-information-in-one-data-frame-1.4.4">Get all selection information in one data frame</a></span></li></ul></li><li><span><a href="#Analyze-immune-selection" data-toc-modified-id="Analyze-immune-selection-1.5">Analyze immune selection</a></span><ul class="toc-item"><li><span><a href="#Examine-all-samples-and-choose-ones-to-retain" data-toc-modified-id="Examine-all-samples-and-choose-ones-to-retain-1.5.1">Examine all samples and choose ones to retain</a></span><ul class="toc-item"><li><span><a href="#Plot-site-level-selection-for-all-samples" data-toc-modified-id="Plot-site-level-selection-for-all-samples-1.5.1.1">Plot site-level selection for all samples</a></span></li><li><span><a href="#Choose-samples-to-retain-based-on-infectivity-remaining" data-toc-modified-id="Choose-samples-to-retain-based-on-infectivity-remaining-1.5.1.2">Choose samples to retain based on infectivity remaining</a></span></li><li><span><a href="#Listing-of-retained-samples" data-toc-modified-id="Listing-of-retained-samples-1.5.1.3">Listing of retained samples</a></span></li></ul></li><li><span><a href="#Compute-serum-average-from-retained-samples" data-toc-modified-id="Compute-serum-average-from-retained-samples-1.5.2">Compute serum average from retained samples</a></span></li><li><span><a href="#Choose-averaging-method-for-downstream-use" data-toc-modified-id="Choose-averaging-method-for-downstream-use-1.5.3">Choose averaging method for downstream use</a></span></li><li><span><a href="#Identify-sites-of-&quot;significant&quot;-selection" data-toc-modified-id="Identify-sites-of-&quot;significant&quot;-selection-1.5.4">Identify sites of "significant" selection</a></span><ul class="toc-item"><li><span><a href="#Cutoff-for-significance" data-toc-modified-id="Cutoff-for-significance-1.5.4.1">Cutoff for significance</a></span></li><li><span><a href="#Identify-significant-sites" data-toc-modified-id="Identify-significant-sites-1.5.4.2">Identify significant sites</a></span></li><li><span><a href="#List-significant-sites-for-each-serum" data-toc-modified-id="List-significant-sites-for-each-serum-1.5.4.3">List significant sites for each serum</a></span></li><li><span><a href="#Get-significant-sites-for-each-serum-group" data-toc-modified-id="Get-significant-sites-for-each-serum-group-1.5.4.4">Get significant sites for each serum group</a></span></li></ul></li><li><span><a href="#Plot-serum-average-selection" data-toc-modified-id="Plot-serum-average-selection-1.5.5">Plot serum-average selection</a></span><ul class="toc-item"><li><span><a href="#Choose-sites-to-zoom-in-on" data-toc-modified-id="Choose-sites-to-zoom-in-on-1.5.5.1">Choose sites to zoom-in on</a></span></li><li><span><a href="#Compact-plots-of-replicate-average-selection" data-toc-modified-id="Compact-plots-of-replicate-average-selection-1.5.5.2">Compact plots of replicate-average selection</a></span></li><li><span><a href="#Compact-plots-showing-each-replicate" data-toc-modified-id="Compact-plots-showing-each-replicate-1.5.5.3">Compact plots showing each replicate</a></span></li><li><span><a href="#Whole-gene-logo-plots-of-replicate-average-selection" data-toc-modified-id="Whole-gene-logo-plots-of-replicate-average-selection-1.5.5.4">Whole-gene logo plots of replicate-average selection</a></span></li></ul></li></ul></li></ul></li><li><span><a href="#Plot-figures-for-paper" data-toc-modified-id="Plot-figures-for-paper-2">Plot figures for paper</a></span></li></ul></div>

# Analysis of serum mutational antigenic profiling
This Python Jupyter notebook analyzes mutational antigenic profiling of serum against virus carrying the A/Perth/2009 (H3N2) HA.

## Configure analysis
We first configure the analysis by importing packages and getting information on the samples that we will analyze.

### Import Python modules / packages
Import modules / packages.
In particular, we use:

 - [plotnine](https://plotnine.readthedocs.io) for ggplot2-like plotting syntax
 - [dmslogo](https://jbloomlab.github.io/dmslogo/) to draw sequence logo plots
 - [dms_tools2](https://jbloomlab.github.io/dms_tools2/) for much of the analysis


```python
import os
import math
import collections
import itertools
import warnings
import subprocess

import yaml
import pandas as pd

from IPython.display import display, HTML

import matplotlib.pyplot as plt
import plotnine as p9  # shorter name instead of import *

import dms_tools2
from dms_tools2.ipython_utils import showPDF
from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY as PALETTE

import dmslogo
```

Turn on interactive matplotlib plotting:


```python
plt.ion()  # turn on interactive plotting
```

Print versions of Bloom lab software:


```python
print(f"Using dms_tools2 version {dms_tools2.__version__}")
print(f"Using dmslogo version {dmslogo.__version__}")
```

Set data frame display options:


```python
pd.set_option('display.max_colwidth', 200)
```

Ignore warnings that can clutter notebook:


```python
warnings.simplefilter('ignore')
```

### Get config variables
Read the variables from the config file.
These variables specify input / output paths and key parameters for analysis:


```python
configfile = 'config.yaml'

with open(configfile) as f:
    config = yaml.safe_load(f)

print(f"Read the following configuration from {configfile}:")
display(HTML(pd.Series(config)
             .to_frame('value')
             .to_html()
             ))
```

### Get sera information
Read information on the sera that are being mapped.

For each serum sample below, we get:
  - *serum*: the abbreviation used for that serum in the experiments
  - *serum_description*: description of the serum
  - *serum_group*: group to which the serum belongs
  - *serum_name*: a short informative name
  - *serum_species*: species from which the serum is derived
  - *serum_vaccination*: any relevant information on vaccination status (some serum are paired pre- and post-vaccination)


```python
with open(config['serum_info']) as f:
    sera = (pd.DataFrame(yaml.safe_load(f))
            .transpose()
            .add_prefix('serum_')
            .rename_axis('serum')
            .reset_index()
            )

assert len(sera) == len(sera['serum'].unique()), 'sera not unique'

print(f"Read the following sera information from {config['serum_info']}:")
display(HTML(sera.to_html(index=False)))
```

### Get list of samples
Read information about all of the samples that we have deep sequenced.

For each sample, we have information on the serum to which it corresponds, the virus library, the date of sequencing, the serum dilution, the percent infectivity, and the location of the FASTQ files on the Hutch server:


```python
samples = pd.read_csv(config['sample_list'])

assert len(samples) == len(samples['sample'].unique()), 'non-unique samples'

print(f"Read the following samples from {config['sample_list']}:")
display(HTML(samples.to_html(index=False)))
```

Check that the serum for all samples are in our set of sera:


```python
unknown_sera = set(samples['serum']) - set(sera['serum'])
if unknown_sera:
    raise ValueError(f"samples include unknown sera: {unknown_sera}")
else:
    print('We have information for all sera used for the samples.')
```

## Download deep sequencing data if needed
The configfile specifies whether we get the data from existing *R1* files on the Hutch server, or download the data from the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA) using [dms_tools2.sra.fastqFromSRA](https://jbloomlab.github.io/dms_tools2/dms_tools2.sra.html):


```python
if config['seq_data_source'] == 'SRA_accession':
    
    # are we using Aspera for rapid downloads?
    if config['ascp'] and config['asperakey']:
        aspera = (config['ascp'], config['asperakey'])
    else:
        aspera = None
        
    # do the downloads
    print(f"Downloading FASTQ files to {config['fastq_dir']} (takes a while)...")
    os.makedirs(config['fastq_dir'], exist_ok=True)
    samples = samples.rename(columns={'sample': 'name',
                                      'SRA_accession': 'run'})
    dms_tools2.sra.fastqFromSRA(
            samples=samples,
            fastq_dump=config['fastq_dump'],
            fastqdir=config['fastq_dir'],
            aspera=aspera,
            ncpus=config['ncpus'],
            )
    samples = samples.rename(columns={'name': 'sample',
                                      'run': 'SRA_accession'})
    print('Completed downloading files.')
        
elif config['seq_data_source'] != 'R1':
    raise ValueError('invalid value of `seq_data_source`')
```

## Align sequencing and count mutations
The samples were sequenced using [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html) to obtain high accuracy.
So we need to process these data to determine the counts of each codon mutation in each sample.

First, create the directory used for the results of this part of the analysis:


```python
os.makedirs(config['countsdir'], exist_ok=True)

print(f"Results from counting mutations go to {config['countsdir']}")
```

### Run `dms2_batch_bcsubamp`
We process the sequencing data by using [dms2_batch_bcsubamp](https://jbloomlab.github.io/dms_tools2/dms2_batch_bcsubamp.html#dms2-batch-bcsubamp) to generate "codon counts" files that give the counts of each codon mutation for each sample.

First, write the batch file used as input by this program:


```python
bcsubamp_batchfile = os.path.join(config['countsdir'], 'batch.csv')

(samples
 .rename(columns={'sample': 'name'})
 [['name', 'R1']]
 .to_csv(bcsubamp_batchfile)
 )

print(f"Creating batch file {bcsubamp_batchfile}")
```

Now run the program:


```python
cmds = ['dms2_batch_bcsubamp',
        '--batchfile', bcsubamp_batchfile,
        '--refseq', config['refseq'],
        '--alignspecs'] + config['alignspecs'].split() + [
        '--R1trim', str(config['R1trim']),
        '--R2trim', str(config['R2trim']),
        '--outdir', config['countsdir'],
        '--fastqdir', config['fastq_dir'],
        '--summaryprefix', 'summary',
        '--ncpus', str(config['ncpus']),
        '--use_existing', config['use_existing'],
        ]

print(f"Running dms2_batch_bcsubamp with this command:\n{' '.join(cmds)}")
subprocess.check_output(cmds)
print('Completed running dms2_batch_bcsubamp.')
```

Confirm that all the expected counts files exist:


```python
assert all(os.path.isfile(f) for f in
           config['countsdir'] + '/' + samples['sample'] + '_codoncounts.csv'
           ), 'missing counts files'
```

### Plot sequencing and mutation counts summaries
Running [dms2_batch_bcsubamp](https://jbloomlab.github.io/dms_tools2/dms2_batch_bcsubamp.html#dms2-batch-bcsubamp) created some summary plots.
The prefix on those plots should be as follows:


```python
countsplotprefix = os.path.join(config['countsdir'], 'summary_')
```

#### Reads and barcodes per sample
Total sequencing reads per sample:


```python
showPDF(countsplotprefix + 'readstats.pdf')
```

Distribution of sequencing reads per barcode on subamplicons:


```python
showPDF(countsplotprefix + 'readsperbc.pdf')
```

Number of barcoded subamplicons that align and have sufficient reads:


```python
showPDF(countsplotprefix + 'bcstats.pdf')
```

#### Coverage across gene
Depth of valid barcoded subamplicons covering each site in the gene:


```python
showPDF(countsplotprefix + 'depth.pdf')
```

#### Mutation frequencies
The average mutation frequency for each sample, stratifying by codon mutation type:


```python
showPDF(countsplotprefix + 'codonmuttypes.pdf')
```

Average mutation frequency per sample, stratifying by number of nucleotide changes per codon mutation:


```python
showPDF(countsplotprefix + 'codonntchanges.pdf')
```

Per-codon mutation frequencies across all sites in gene for each sample:


```python
showPDF(countsplotprefix + 'mutfreq.pdf')
```

#### Check for oxidative damage
Sometimes there is oxidative damage which manifests as an enrichment of `G`->`T` and `C`->`A` mutations among the single-nucleotide codon mutations.
Check for this by plotting frequencies of different single-nucleotide mutation types:


```python
showPDF(countsplotprefix + 'singlentchanges.pdf')
```

### Renumber to standard HA numbering
The above alignments use sequential 1, 2, ... numbering of the codons.
This is not the standard numbering scheme used for HA, so we use [dms_tools2.utils.renumbSites](https://jbloomlab.github.io/dms_tools2/dms_tools2.utils.html#dms_tools2.utils.renumberSites) to renumber to the standard HA numbering scheme:


```python
dms_tools2.utils.renumberSites(
    renumbfile=config['renumbering_scheme'],
    infiles=list(config['countsdir'] + '/' + samples['sample'] +
                 '_codoncounts.csv'),
    missing='drop',
    outdir=config['renumbcountsdir'])

assert all(os.path.isfile(f) for f in
           config['renumbcountsdir'] + '/' + samples['sample'] +
           '_codoncounts.csv'
           ), 'missing renumbered counts files'

print(f"Renumbered codon counts are in {config['renumbcountsdir']}")
```

## Compute immune selection on mutations
We will now quantify the immune selection on each mutation by comparing its frequency in each serum-selected sample to an appropriate mock-selected control.

### Samples to compare for each selection
In order to quantify the immune selection, we need to compare each selected sample to the appropriate controls.
Specifically, for each selection, we define three samples:
 - *sel*: the immune-selected sample
 - *mock*: the appropriate mock-selected control for that date and library
 - *err*: the appropriate wildtype plasmid control to estimate sequencing error rates.
 
We call each combination of three such samples a "selection".
For each selection, we also record:
  - *serum_name_formatted*: a version of *serum_name* that also indicates species if non-human
  - *name*: a string indicating the library and the percent infectivity remaining. So for instance, *lib1-5* indicates library 1 with 5% infectivity remaining.
  
Below we construct the data frame with this information on the selections:


```python
selections = (

    # get all immune selected (sel) samples
    samples
    .query('(serum != "mock") & (serum != "plasmid")')
    .rename(columns={'sample': 'sel'})

    # add mock sample for that date and library
    .merge(samples
           .query('serum == "mock"')
           [['sample', 'library', 'date']]
           .rename(columns={'sample': 'mock'})
           )

    # add plasmid sample as error control (err) for that date and library
    .merge(samples
           .query('serum == "plasmid"')
           [['sample', 'library', 'date']]
           .rename(columns={'sample': 'err'})
           )

    # add information about sera
    .merge(sera, validate='many_to_one')

    # add columns with frac survive and informative names for serum and samples
    .assign(
        libfracsurvive=lambda x: x['percent_infectivity'] / 100,
        serum_name_formatted=lambda x:
            x['serum_species'].map(lambda s: '' if pd.isnull(s) or s == 'human'
                                   else s + '-') + x['serum_name'],
        name_formatted=lambda x:
            x['library'] + ', ' + x['percent_infectivity'].apply(
                dms_tools2.utils.sigFigStr, nsig=2) + '% infectivity',
        name=lambda x:
            x['library'] + '-' + x['percent_infectivity'].apply(
                dms_tools2.utils.sigFigStr, nsig=2)
        )

    # drop unneeded columns
    .drop(['R1'], axis='columns')

    # re-order columns a bit so key ones are displayed first
    .set_index(['serum_name_formatted', 'name', 'sel', 'mock', 'err',
                'libfracsurvive'])
    .reset_index()
    )

# make sure no duplicated serum / names
assert len(selections) == len(selections.groupby(['serum_name_formatted',
                                                  'name']))

print(f"Tabulated information for {len(selections)} selections:")
display(HTML(selections.to_html(index=False)))
```

### Choose measure of immune selection
In our prior work, we have used two different measures of the selection on each mutation:
 - *differential selection (diffsel)*: essentially the log enrichment of each mutation relative to wildtype in the immune-selected versus mock sample; [see here for a detailed description](https://jbloomlab.github.io/dms_tools2/diffsel.html).
 - *fraction surviving above average (fracsurvive)*: the estimated fraction of virions with each mutation that survive immune selection **above the library average**; [see here for a detailed description](https://jbloomlab.github.io/dms_tools2/fracsurvive.html).
 
There are programs in [dms_tools2](https://jbloomlab.github.io/dms_tools2) to compute each of these measures of selection, namely [dms_batch_diffsel](https://jbloomlab.github.io/dms_tools2/dms2_batch_diffsel.html) or [dms2_batch_fracsurvive](https://jbloomlab.github.io/dms_tools2/dms2_batch_fracsurvive.html).

For this notebook, we will just choose **one** of these measures to compute.
The reason is that we will use just one measure in the paper, and the notebook is less cluttered if we show just one measure of selection.
Note however that the notebook is set up so that you can simply switch the value of the `seltype` variable defined immediately below and then run all the results for the other measure of selection if you want to compare things:


```python
seltype = 'diffsel'  # use this measure of immune selection
```

We also specify the extra arguments that must be passed to the [dms_batch_diffsel](https://jbloomlab.github.io/dms_tools2/dms2_batch_diffsel.html) or [dms2_batch_fracsurvive](https://jbloomlab.github.io/dms_tools2/dms2_batch_fracsurvive.html) program for each selection metric:


```python
# extra arguments when computing selection
addargs = {'diffsel': [],
           'fracsurvive': ['--aboveavg', 'yes']  # fracsurvive above average
           }
```

The selection is computed at the mutation level.
But for many plots, we want to summarize the selection at the site level.
So we need to define what site-level metric we use:


```python
site_metric = {'diffsel': 'positive_diffsel',
               'fracsurvive': 'avgfracsurvive'
               }
```

### Compute immune selection
Now we actually run the appropriate program to compute the immune selection.
We then add to our `selections` data frame the name of the files holding the computed site (*site*) and mutation (*mut*) level selection for each sample.
  
The next cell does all of this:


```python
prog = f"dms2_batch_{seltype}"
outdir = config[f"{seltype}dir"]
os.makedirs(outdir, exist_ok=True)

# write batch file used by program
batchfile = os.path.join(outdir, 'batch.csv')
(selections
 .rename(columns={'serum_name_formatted': 'group'})
 .to_csv(batchfile, index=False)
 )

cmds = [prog,
        '--summaryprefix', 'summary',
        '--batchfile', batchfile,
        '--outdir', outdir,
        '--indir', config['renumbcountsdir'],
        '--use_existing', config['use_existing'],
        '--ncpus', str(config['ncpus'])
        ] + addargs[seltype]

print(f"Computing {seltype} using {prog} with command:\n{' '.join(cmds)}")
subprocess.check_output(cmds)

selfilecols = []
for selfile in ['mut' + seltype, 'site' + seltype]:
    selfilecol = selfile + '_file'
    selfilecols.append(selfilecol)
    selections[selfilecol] = (outdir + '/' + selections['serum_name_formatted']
                              + '-' + selections['name'] + '_' +
                              selfile + '.csv')
    assert all(selections[selfilecol].map(os.path.isfile)), 'missing files'
    print(f"Created {len(selections[selfilecol])} {selfile} files, adding to "
          f"`selections` data frame in column {selfilecol}")
```

### Get all selection information in one data frame
For further processing, we want to create a dataframe that holds all of the selection information at the site and mutation levels for all samples.
We create such a dataframe, *sel_df*, by reading the files in *selections* into the data frame using [dms_tools2.diffsel.df_read_filecols](https://jbloomlab.github.io/dms_tools2/dms_tools2.diffsel.html#dms_tools2.diffsel.df_read_filecols):


```python
sel_df = dms_tools2.diffsel.df_read_filecols(selections, selfilecols)
```

Now *sel_df* is a very large data frame, but it has all the information we want to plot.
Here are the first few rows:


```python
print(f"sel_df has {len(sel_df)} rows. Here are the first few:")
display(HTML(sel_df.head(n=5).to_html(index=False)))
```

## Analyze immune selection

### Examine all samples and choose ones to retain
For some serum, we have multiple samples.
These include the samples from multiple library replicates, but sometimes also several samples for the same library replicate at different serum concentrations that yield different percent viral infectivity remaining.

Which of these samples do we want to retain?
We probably want just one sample per library (so we don't over-weight some libraries in our results), and we need to select the serum concentration (percent viral infectivity remaining) that is "best" in terms of giving reproducible results between replicates with maximal signal to noise.

#### Plot site-level selection for all samples
We first plot the site-level selection for all samples.
To do this, we look over all sites and use the `facet_plot` command of [dmslogo](https://jbloomlab.github.io/dmslogo/) to plot the site-level selection for all samples for each serum.
We also compute the correlation between samples for each serum:


```python
for serum_name, serum_sel_df in sel_df.groupby('serum_name_formatted'):

    print(f"\n\n******************* {serum_name} *******************")

    fig, ax = dmslogo.facet_plot(
                serum_sel_df,
                x_col='isite',
                gridcol_col='name_formatted',
                show_col=None,
                wspace=0.6,
                draw_line_kwargs=dict(
                        xtick_col='site',
                        height_col=site_metric[seltype],
                        ylabel=site_metric[seltype],
                        )
                )
    display(fig)
    plt.close(fig)

    corr_df = (serum_sel_df
               .rename(columns={'name_formatted': 'sample'})
               .pivot_table(values=site_metric[seltype],
                            columns='sample',
                            index=['site'])
               .corr()
               .round(3)
               )
    display(HTML(corr_df.to_html()))
```

#### Choose samples to retain based on infectivity remaining
From the plots and correlations immediately above, several things are apparent:
 1. Different samples with the **same** virus library tend to be substantially more correlated than samples between different virus libraries.
    This makes sense, as some of the noise probably comes from library composition differences in terms of variant abundance.
    This fact implies that we want to keep just one (or the same number) of samples from each library, so that when we take means / medians across samples they are not disproportionately influenced by a single library.
 2. The "best" concentration of serum appears to be one that leaves about about 1% to 5% of the viral infectivity remaining. 
    This is judged by the fact that samples with this percent infectivity remaining tend to look "cleanest" in the line plots of site differential selection, and also tend to correlate best with samples from another library.
    
Based on the above, we will filter our samples down to those that have closest to the target amount of viral infectivity remaining.
Because we qualitatively estimated above that the best results are obtained when between 1% to 5% of infectivity is remaining, we will target 2% infectivity.
Then for any library / serum with multiple samples, we will retain the one that is closest (in log space) to that target infectivity.
We do this in log space because in linear space, 0% infectivity would be closer to 2% infectivity as would 6% infectivity--even though the latter is clearly preferable if our target is 2%.


```python
target_infectivity = 2
print(f"Choosing samples closest to {target_infectivity:.2f}% infectivity.")
```

Mark the samples to retain by adding a *retained* column to both the *selections* and *sel_df* data frames, and flagging it as *True* only for the sample for that serum and library that has the smallest log-space distance to the target infectivity:


```python
selections = (
    selections
    .assign(retained=lambda x: x.assign(dist=lambda y: (y.percent_infectivity /
                                                        target_infectivity)
                                        .apply(math.log).apply(abs))
            .groupby(['serum_name_formatted', 'library'])
            ['dist']
            .transform(lambda y: y <= y.min())
            )
    )
print(f"Retaining {len(selections.query('retained'))} of {len(selections)}")

sel_df = sel_df.merge(selections[['serum_name_formatted', 'name', 'retained']],
                      validate='many_to_one')
```

Plot the samples to retain and their percent infectivity.
In the plot below, the dashed horizontal line indicates the target percent infectivity, and the colors of the points indicate which sample we retained for each serum / library:


```python
p = (p9.ggplot(selections.assign(xlabel=lambda x: (x['serum_name_formatted'] +
                                                   ', ' + x['library'])),
               p9.aes('xlabel', 'percent_infectivity', color='retained')) +
     p9.geom_point(size=3, alpha=0.7) +
     p9.theme(
         axis_text_x=p9.element_text(angle=90),
         figure_size=(0.25 * len(selections.groupby(['serum_name_formatted',
                                                     'library'])), 2.5)
         ) +
     p9.scale_y_log10(name='percent infectivity') +
     p9.xlab('') +
     p9.scale_color_manual(values=PALETTE) +
     p9.geom_hline(yintercept=target_infectivity, linetype='dashed',
                   alpha=0.7, color=PALETTE[2])
     )

_ = p.draw()
```

#### Listing of retained samples
Here is a small table listing the samples that we retained for each serum, and their percent infectivity remaining:


```python
display(HTML(selections
             .assign(serum_name=lambda x:
                     pd.Categorical(x['serum_name'],
                                    sera['serum_name'],
                                    ordered=True))
             .query('retained')
             .sort_values('serum_name')
             [['serum_group', 'serum_name_formatted', 'library',
               'percent_infectivity']]
             .pivot_table(index=['serum_group', 'serum_name_formatted'],
                          columns='library',
                          values='percent_infectivity')
             .to_html()
             ))
```

### Compute serum average from retained samples
We now compute the "average" selection among the retained samples for each serum.

First, confirm that we have retained just one sample per serum / library:


```python
assert all(len(group) == 1 for _, group in
           (selections
            .query('retained')
            .groupby(['serum_name_formatted', 'library'])
            )
           )
```

We will compute two types of "average" selection, the mean and the median across the samples:


```python
avg_types = ['mean', 'median']
```

Do a sanity check and make sure none of our libraries are named one of the average types:


```python
assert not set(avg_types).intersection(set(selections['library']))
assert not set(avg_types).intersection(set(selections['name_formatted']))
```

Now loop over all sera and compute the average selection for all retained samples for that sera.
Note that the averages are computed on the mutation-level data, and then the site data are computed from those averaged mutation-level data.
The averages (along with the samples used to compute these averages) are then added to a new data frame similar to *selections* that is called *avg_selections*:


```python
# function to compute mutation-level average
mut_avg_func = {'diffsel': dms_tools2.diffsel.avgMutDiffSel,
                'fracsurvive': dms_tools2.fracsurvive.avgMutFracSurvive}

# function to convert mutation-level values to site-level values
mut_to_site_func = {'diffsel': dms_tools2.diffsel.mutToSiteDiffSel,
                    'fracsurvive': dms_tools2.fracsurvive.mutToSiteFracSurvive}

avg_selections = []
for serum_name_formatted, group in (
            selections
            .query('retained')
            [['serum_name_formatted', 'library', 'name_formatted'] +
             list(sera.columns) + selfilecols]
            .groupby('serum_name_formatted')
            ):

    avg_selections.append(group)

    for avg_type in avg_types:
        # build row of data frame with average
        avg_row = group.iloc[0].to_dict(into=collections.OrderedDict)
        avg_row['library'] = avg_type
        avg_row['name_formatted'] = avg_type

        avgdir = config[f"avg{seltype}dir"]
        os.makedirs(avgdir, exist_ok=True)

        avg_row[f"mut{seltype}_file"] = (f"{avgdir}/{serum_name_formatted}-"
                                         f"mut{seltype}-{avg_type}.csv")
        (mut_avg_func[seltype](group[f"mut{seltype}_file"], avg_type)
         .to_csv(avg_row[f"mut{seltype}_file"], index=False))

        avg_row[f"site{seltype}_file"] = (f"{avgdir}/{serum_name_formatted}-"
                                          f"site{seltype}-{avg_type}.csv")
        (mut_to_site_func[seltype](pd.read_csv(avg_row[f"mut{seltype}_file"]))
         .to_csv(avg_row[f"site{seltype}_file"], index=False))

        avg_row = pd.Series(avg_row).to_frame().transpose()
        assert all(avg_row.columns == group.columns)
        avg_selections.append(avg_row)

# put avg_selections in data frame, sort to match config['sera']
avg_selections = (pd.concat(avg_selections)
                  .assign(serum_name=lambda x:
                          pd.Categorical(x['serum_name'],
                                         sera['serum_name'],
                                         ordered=True))
                  .sort_values(['serum_name', 'library'])
                  .assign(serum_name_formatted=lambda x:
                          pd.Categorical(x['serum_name_formatted'],
                                         x['serum_name_formatted'].unique(),
                                         ordered=True))
                  .reset_index(drop=True)
                  )
```

Now the `avg_selections` data frame lists the files giving all the retained library replicates plus the means and medians calculated from them these replicates.

Now we create the data frame `avg_sel_df` which actually holds the site- and mutation-level averages for all sera as well as the samples (one replicate per library) that we used to compute these averages:


```python
avg_sel_df = (dms_tools2.diffsel.df_read_filecols(avg_selections, selfilecols)
              # preserve order of sera as in `avg_selections`
              .assign(serum_name_formatted=lambda x:
                      pd.Categorical(x['serum_name_formatted'],
                                     (avg_selections['serum_name_formatted']
                                      .unique()),
                                     ordered=True))
              )
```

This `avg_sel_df` data frame differs from `sel_df` only in that it includes the averages as a library type, and only has the retained replicates for each library.

### Choose averaging method for downstream use
For the summary plots below, we need to choose whether to represent our "averages" using the mean or the median.
This is done in the next cell.
We are using the *median* as our measure of the average; to instead use the *mean* simply switch the cell below to define `avg_type` as *mean* rather than *median*.

**(This choice of median over mean should be re-visited once we have replicates for all three libraries, but for 1 and 2 libraries the mean and the median are the same so it doesn't matter).**


```python
avg_type = 'median'
unused_avg_types = [x for x in avg_types if x != avg_type]
```

### Identify sites of "significant" selection
We want to identify the sites that are under "significant" immune selection.
The reason is that we can then zoom in on these sites in logo plots.

In [dms_tools2.plot.findSigSel](https://jbloomlab.github.io/dms_tools2/dms_tools2.plot.html#dms_tools2.plot.findSigSel) function, we have defined a heuristic way to do this.
Essentially, this function uses robust regression to fit a gamma distribution to the selection values for each site, and then identifies those that are "outliers" of high selection.

#### Cutoff for significance
First, we define a cutoff for what constitutes significant:


```python
fdr_cutoff = 0.05
```

#### Identify significant sites
Now we use [dms_tools2.plot.findSigSel](https://jbloomlab.github.io/dms_tools2/dms_tools2.plot.html#dms_tools2.plot.findSigSel) to get a dataframe (`sigsites_df`) that lists the "significant" sites for each serum.
Note that the cell below also saves plots showing the fit gamma distribution (you can inspect these separately if you want to look in more detail):


```python
plotfile_template = os.path.join(config[f"avg{seltype}dir"],
                                 'sigsites_{serum}.pdf')

print(f"Identifying sites of significant selection at a FDR of {fdr_cutoff}.\n"
      f"Plots of distribution fitting saved as {plotfile_template}")

sigsites_df = []
for serum_name_formatted, group in (
        avg_sel_df
        .query('library == @avg_type')
        [['serum_group', 'serum_name_formatted', 'isite', 'site',
          site_metric[seltype]]]
        .drop_duplicates()
        .groupby('serum_name_formatted')
        ):
    plotfile = plotfile_template.format(serum=serum_name_formatted)
    df, cutoff, gamma_params = dms_tools2.plot.findSigSel(
            group,
            site_metric[seltype],
            plotfile,
            fdr=fdr_cutoff,
            title=serum_name_formatted
            )
    sigsites_df.append(df)

sigsites_df = pd.concat(sigsites_df, ignore_index=True).query('sig')

print('Here are the first few rows of sigsites_df:')
display(HTML(sigsites_df.head(n=4).to_html(index=False)))
```

#### List significant sites for each serum
Now display lists of the significant sites for each serum:


```python
display(HTML(sigsites_df
             .sort_values('isite')
             .assign(nsites=1)
             .groupby('serum_name_formatted')
             .aggregate({'site': lambda x: ', '.join(list(x)),
                         'nsites': 'sum'})
             .rename(columns={'site': 'significant sites',
                              'nsites': 'number of sites'})
             .to_html()
             ))
```

#### Get significant sites for each serum group
In the analyses below, we may want to plot the sites that are significant within at least one serum sample for each serum group.
Therefore, we determine the significant sites within each serum group:


```python
sigsites_by_serumgroup = (
    sigsites_df
    [['serum_group', 'isite', 'site']]
    .drop_duplicates()
    .sort_values('isite')
    .groupby('serum_group')
    .aggregate(lambda x: list(x))
    .assign(nsites=lambda x: x['site'].apply(len))
    .sort_values('nsites')
    )

display(HTML(sigsites_by_serumgroup.to_html()))
```

### Plot serum-average selection
Now we plot the average (across libraries) selection for each serum.

#### Choose sites to zoom-in on
In the plots, we will zoom in on important sites using logo plots.
The reason that we zoom in on just a subset of sites is to keep the logo plots relatively compact.

The sites we zoom in on will be those identified above as being under "significant" selection.
There are two ways we can do the zooming:
 1. We can zoom in on separate sites for each serum group.
    In this case, for each group we only zoom on sites of significant selection for at least one serum in that group.
    This might be preferable if each serum group targets very different sites.
    To do this, set *zoom_combine_serum_groups* to *False*.
 2. We can zoom in on the same sites for all serum groups.
    In this case, for each group we zoom on all sites that are significant for any serum in any group.
    This might be preferable if we want to compare across sites.
    To do this, set *zoom_combine_serum_groups* to *True*.


```python
zoom_combine_serum_groups = False
```

We also may want to "pad" the sites that we zoom in on by zooming in on this many sites before and after each significant site.
This is useful if you want to give some context to the zoomed sites.
Below, set the amount of padding you want around each significant site; a value of 0 means no padding:


```python
zoom_pad = 0
```

Now build a dict, *zoom_sites* that is keyed first by *serum_group* and then by *isite* and *site*, with values being the list of sites to zoom for that serum group:


```python
zoom_sites = {}
for tup in sigsites_by_serumgroup.reset_index().itertuples(index=False):
    if zoom_combine_serum_groups:
        isites = set(itertools.chain.from_iterable(
                     sigsites_by_serumgroup['isite']))
    else:
        isites = set(tup.isite)
    for isite in list(isites):
        for pad in range(-zoom_pad, zoom_pad + 1):
            isites.add(isite + pad)
    # only keep valid sites
    isites = isites.intersection(set(avg_sel_df['isite']))
    isites = sorted(isites)  # get as sorted list
    # get sites corresponding to each isite
    sites = (avg_sel_df
             [['isite', 'site']]
             .drop_duplicates()
             .sort_values('isite')
             .query('isite in @isites')
             ['site']
             .tolist()
             )
    zoom_sites[tup.serum_group] = {'isite': isites,
                                   'site': sites,
                                   'nsites': len(sites)}
```

Here are the sites we will zoom in on for each serum group:


```python
display(HTML(pd.DataFrame.from_dict(zoom_sites, orient='index').to_html()))
```

Add a column (`zoom_site`) to `avg_sel_df` that indicates which sites to zoom in on:


```python
avg_sel_df = pd.concat(
        [df.assign(zoom_site=lambda x: x['isite'].isin(
                            zoom_sites[serum_group]['isite']))
         for serum_group, df in avg_sel_df.groupby('serum_group')],
        ignore_index=True
        )
```

#### Compact plots of replicate-average selection
For each group of sera we make line plots that show the site-level selection and logo plots that zoom in on mutations at the sites of significant selection.
We make these plots using the `facet_plot` command of [dmslogo](https://jbloomlab.github.io/dmslogo/).

We want to label the logo plots with site numbers **and** wildtype residue, so first we create a column that contains this information:


```python
avg_sel_df = avg_sel_df.assign(site_label=lambda x: x['wildtype'] +
                               x['site'].astype('str'))
```

Now we make the line and logo plots.
We also save PDF versions of each plot:


```python
for serum_group, df in avg_sel_df.groupby('serum_group'):

    plotfile = os.path.join(config[f"avg{seltype}dir"],
                            f"{serum_group}_avg.pdf")
    print(f"\n\n{'*' * 72}\nSerum group {serum_group}, saving to {plotfile}\n")

    fig, axes = dmslogo.facet_plot(
            data=df.query('library == @avg_type'),
            x_col='isite',
            show_col='zoom_site',
            gridrow_col='serum_name_formatted',
            share_xlabel=True,
            share_ylabel=True,
            wspace=0.6,
            draw_line_kwargs=dict(
                    height_col=site_metric[seltype],
                    xtick_col='site',
                    ylabel=f"immune selection ({seltype})",
                    ),
            draw_logo_kwargs=dict(
                    letter_col='mutation',
                    letter_height_col=f"mut{seltype}",
                    xtick_col='site_label',
                    xlabel='site',
                    ylabel=f"immune selection ({seltype})",
                    clip_negative_heights=True,
                    ),
            )
    display(fig)
    fig.savefig(plotfile)
    plt.close(fig)
```

#### Compact plots showing each replicate
The above plots show the average of the retained replicates.
It can also be helpful to look at each replicate that went into the average separately.

The next cell makes compact plots for each replicate separately:


```python
for serum_group, df in avg_sel_df.groupby('serum_group'):

    plotfile = os.path.join(config[f"avg{seltype}dir"],
                            f"{serum_group}_reps.pdf")
    print(f"\n\n{'*' * 72}\nSerum group {serum_group}, saving to {plotfile}\n")

    fig, axes = dmslogo.facet_plot(
            data=df.query('library not in @avg_types'),
            x_col='isite',
            show_col='zoom_site',
            gridrow_col='serum_name_formatted',
            gridcol_col='library',
            share_xlabel=True,
            share_ylabel=True,
            wspace=0.6,
            draw_line_kwargs=dict(
                    height_col=site_metric[seltype],
                    xtick_col='site',
                    ylabel=f"immune selection ({seltype})",
                    ),
            draw_logo_kwargs=dict(
                    letter_col='mutation',
                    letter_height_col=f"mut{seltype}",
                    xtick_col='site_label',
                    xlabel='site',
                    ylabel=f"immune selection ({seltype})",
                    clip_negative_heights=True,
                    ),
            )
    display(fig)
    fig.savefig(plotfile)
    plt.close(fig)
```

#### Whole-gene logo plots of replicate-average selection
Finally, we make whole-gene logo plots for each serum that shows the replicate-average selection for **all** sites.
We make these whole-gene plots using [dms2_logoplot](https://jbloomlab.github.io/dms_tools2/dms2_logoplot.html).
They are too large to be useful to show visually in this notebook, but the cell below gives the name of the PDF holding each logo plot so you can examine them individually if you'd like.


```python
outdir = config[f"avg{seltype}dir"]  # save plots here

for tup in (avg_selections
            .query('library == @avg_type')
            .itertuples(index=False)
            ):
    name = getattr(tup, "serum_name_formatted")
    plotfile = os.path.join(outdir, f"{name}_{seltype}.pdf")
    datafile = getattr(tup, f"mut{seltype}_file")
    cmds = ['dms2_logoplot',
            '--outdir', outdir,
            f"--{seltype}", datafile,
            '--name', name,
            '--nperline', '71',
            '--overlay1', datafile, 'wildtype', 'wildtype',
            '--underlay', 'yes', 
            '--restrictdiffsel', 'positive',
            '--use_existing', config['use_existing'],
            ]
    print(f"Plotting {name} to {plotfile}")
    subprocess.check_output(cmds)
    assert os.path.isfile(plotfile)
```

# Plot figures for paper

In the following section, we will generate the plots and figures that will go into the paper.

First, we will plot the profiles for the serum with monoclonal Ab spike-in as well as the profiles for the individual serum and mAb. Significant sites of escape for the mAb will be in one color, and sites for the serum will be in another color.


```python
spike_samples = ['antibody-4F03', '2009-age-65']
spike_sample_isites = {}

for s in spike_samples:
    l = list((sigsites_df.query('serum_name_formatted == @s')
              [['serum_name_formatted', 'isite', 'site']]
              .drop_duplicates()
              .sort_values('isite'))['isite'])
    spike_sample_isites[s] = l
```


```python
spike_zoom_sites = [x for v in spike_sample_isites.values() 
                    for x in v]
```


```python
avg_sel_df = pd.concat(
        [df.assign(spike_zoom_site=lambda x: x['isite'].isin(
                            spike_zoom_sites))
         for serum_group, df in avg_sel_df.groupby('serum_group')],
        ignore_index=True
        )
```


```python
avg_sel_df['color'] = (avg_sel_df['isite']
                       .apply(lambda x: PALETTE[2] 
                              if x in spike_sample_isites['antibody-4F03'] 
                              else (PALETTE[1] if x in spike_sample_isites['2009-age-65'] 
                                    else 'gray')))
```


```python
spike_ids = {'antibody-4F03': 'mAb', 
             '2009-age-65': 'serum', 
             '2009-age-65-with-low-4F03': 'serum with low mAb', 
             '2009-age-65-with-mid-4F03': 'serum with intermediate mAb', 
             '2009-age-65-with-hi-4F03': 'serum with high mAb'}

avg_sel_df['spike_id'] = (avg_sel_df['serum_name_formatted']
                            .apply(lambda x: spike_ids[x] 
                                   if x in spike_ids.keys() 
                                   else 'none'))

avg_sel_df['spike_group'] = (avg_sel_df['serum_name_formatted']
                                .apply(lambda x: 'spike_group' 
                                   if x in spike_ids.keys() 
                                   else 'not spike_group'))
```


```python
plotfile = os.path.join(config[f"avg{seltype}dir"], "spikein_avg.pdf")

fig, axes = dmslogo.facet_plot(
        data=avg_sel_df.query('spike_group == "spike_group" and library == @avg_type'), 
        x_col='isite', 
        show_col='spike_zoom_site', 
        gridrow_col='spike_id', 
        share_xlabel=True, 
        share_ylabel=True, 
        wspace=0.6,  
        draw_line_kwargs=dict(
                height_col=site_metric[seltype], 
                xtick_col='site', 
                ylabel=f"immune selection ({seltype})",
                ), 
        draw_logo_kwargs=dict(
                letter_col='mutation', 
                letter_height_col=f"mut{seltype}", 
                color_col = 'color',
                xtick_col='site_label', 
                xlabel='site', 
                ylabel=f"immune selection ({seltype})", 
                clip_negative_heights=True,
                ),
        )

display(fig)
fig.savefig(plotfile)
plt.close(fig)
```


```python

```
