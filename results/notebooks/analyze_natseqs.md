
## Setup notebook
Import Python modules:


```python
import collections
import gzip
import multiprocessing
import os
import re
import subprocess
import warnings

import Bio.SeqIO
import Bio.SeqRecord

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import yaml

import dms_tools2
from dms_tools2 import AAS
from dms_tools2.ipython_utils import showPDF
from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY as PALETTE
print(f"Using dms_tools2 version {dms_tools2.__version__}")
```

    Using dms_tools2 version 2.4.12


Filter warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read configuration:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Get output directory and create if needed:


```python
resultsdir = config['natseqs_dir']
os.makedirs(resultsdir, exist_ok=True)
```

CPUs to use:


```python
ncpus = min(config['ncpus'], multiprocessing.cpu_count())
```

## Align and filter sequences
Read all of the natural sequences and add the reference sequence (the HA used in deep mutational scanning):


```python
print(f"Processing sequences from {config['natseqs']}")
with gzip.open(config['natseqs'], 'rt') as f:
    seqs = list(Bio.SeqIO.parse(f, 'fasta'))
print(f"There are {len(seqs)} sequences.")

print(f"Adding reference sequence in {config['refseq']}")
refseq = Bio.SeqIO.read(config['refseq'], 'fasta')
seqs.append(refseq)
```

    Processing sequences from data/human_H3N2_HA_2007-2018.fasta.gz
    There are 19094 sequences.
    Adding reference sequence in data/Perth09_HA_reference.fa


Use [phydms_prepalignment](http://jbloomlab.github.io/phydms/phydms_prepalignment.html) to align the coding sequences from the proteins, only keeping sequences that are unique at the protein level.
We also only keep sequences with $\ge$92% identity with the reference sequence, as this appears to be a good cutoff between true human H3N2 from this timeframe and mis-annotated sequences or ones from non-human H3 lineages:


```python
print('Using following version of phydms:\n' +
      subprocess.check_output(['phydms_prepalignment', '--version']).decode('utf-8'))

inseqsfile = os.path.join(resultsdir, 'all_seqs.fasta')
Bio.SeqIO.write(seqs, inseqsfile, 'fasta')

alignmentfile = os.path.join(resultsdir, 'all_alignment.fasta')

if config['use_existing'] == 'yes' and os.path.isfile(alignmentfile):
    print(f"Using existing alignment in {alignmentfile}")
else:
    print(f"Building alignment and writing to {alignmentfile}...")
    _ = subprocess.check_call([
        'phydms_prepalignment',
        inseqsfile,
        alignmentfile,
        refseq.name,
        '--minuniqueness', '1',
        '--minidentity', '0.92',
        ])
    print("Bulding of alignment complete.")
```

    Using following version of phydms:
    phydms_prepalignment 2.3.1
    
    Using existing alignment in results/natseqs/all_alignment.fasta


Look at plot showing divergence from "reference" seqeuence:


```python
showPDF(os.path.splitext(alignmentfile)[0] + '.pdf')
```


![png](analyze_natseqs_files/analyze_natseqs_15_0.png)


Read alignment as proteins:


```python
alignment = [Bio.SeqRecord.SeqRecord(seq.seq.translate(gap='-'), name=seq.name)
             for seq in Bio.SeqIO.parse(alignmentfile, 'fasta')
             if seq.name != refseq.name]
refprot = refseq.seq.translate(to_stop=True)
assert all(len(refprot) == len(seq) for seq in alignment)

print(f"Read all {len(alignment)} aligned proteins.")
```

    Read all 5382 aligned proteins.


## Get amino-acid frequencies and other relevant information
Get amino-acid frequencies at each site for each year, breaking into 6-month intervals.

First, we partition sequences into these 6-month intervals:


```python
unmatched = 0
alignment_by_year = collections.defaultdict(list)
for seq in alignment:
    datematch = re.search('_(?P<year>\d{4})/(?P<month>\d{1,2})/\d*_HA$', seq.name)
    if not datematch:
        unmatched += 1
    else:
        year = int(datematch.group('year'))
        month = int(datematch.group('month'))
        if 1 <= month <= 3:
            year = year
        elif 3 < month <= 9:
            year = year + 0.5
        elif 9 < month <= 12:
            year = year + 1
        else:
            raise ValueError('invalid month')
        alignment_by_year[year].append(str(seq.seq))

print('Here are number of sequences in each half-year:')
display(HTML(pd.DataFrame.from_records([(year, len(seqs)) for 
                                        year, seqs in alignment_by_year.items()],
                                       columns=['year', 'nsequences'])
             .sort_values('year')
             .to_html(index=False)
             ))
```

    Here are number of sequences in each half-year:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>year</th>
      <th>nsequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2007.0</td>
      <td>85</td>
    </tr>
    <tr>
      <td>2007.5</td>
      <td>67</td>
    </tr>
    <tr>
      <td>2008.0</td>
      <td>111</td>
    </tr>
    <tr>
      <td>2008.5</td>
      <td>15</td>
    </tr>
    <tr>
      <td>2009.0</td>
      <td>76</td>
    </tr>
    <tr>
      <td>2009.5</td>
      <td>158</td>
    </tr>
    <tr>
      <td>2010.0</td>
      <td>27</td>
    </tr>
    <tr>
      <td>2010.5</td>
      <td>105</td>
    </tr>
    <tr>
      <td>2011.0</td>
      <td>260</td>
    </tr>
    <tr>
      <td>2011.5</td>
      <td>98</td>
    </tr>
    <tr>
      <td>2012.0</td>
      <td>303</td>
    </tr>
    <tr>
      <td>2012.5</td>
      <td>159</td>
    </tr>
    <tr>
      <td>2013.0</td>
      <td>500</td>
    </tr>
    <tr>
      <td>2013.5</td>
      <td>76</td>
    </tr>
    <tr>
      <td>2014.0</td>
      <td>161</td>
    </tr>
    <tr>
      <td>2014.5</td>
      <td>104</td>
    </tr>
    <tr>
      <td>2015.0</td>
      <td>712</td>
    </tr>
    <tr>
      <td>2015.5</td>
      <td>208</td>
    </tr>
    <tr>
      <td>2016.0</td>
      <td>274</td>
    </tr>
    <tr>
      <td>2016.5</td>
      <td>199</td>
    </tr>
    <tr>
      <td>2017.0</td>
      <td>731</td>
    </tr>
    <tr>
      <td>2017.5</td>
      <td>187</td>
    </tr>
    <tr>
      <td>2018.0</td>
      <td>490</td>
    </tr>
    <tr>
      <td>2018.5</td>
      <td>82</td>
    </tr>
    <tr>
      <td>2019.0</td>
      <td>90</td>
    </tr>
  </tbody>
</table>


Build up data frame giving amino-acid counts at each site in each year, then compute frequencies and add site numbers in H3 numbering:


```python
aafreqs_df = []
for year, year_alignment in alignment_by_year.items():
    for r, wt in enumerate(refprot):
        counts = collections.Counter([seq[r] for seq in year_alignment])
        aafreqs_df.append(
            pd.DataFrame({
                'year': year,
                'isite': r,
                'wildtype': wt,
                'mutation': [tup[0] for tup in sorted(counts.items())],
                'natural_count': [tup[1] for tup in sorted(counts.items())],
                })
            .assign(natural_frequency=lambda x: x['natural_count'] / len(year_alignment)) 
            )

aafreqs_df = pd.concat(aafreqs_df)

print(f"Adding H3 site numbering from {config['renumbering_scheme']}")

aafreqs_df = aafreqs_df.merge(
                pd.read_csv(config['renumbering_scheme'])
                .rename(columns={'original': 'isite', 'new': 'site'})
                .assign(isite=lambda x: x['isite'] - 1)
                )
```

    Adding H3 site numbering from data/H3renumbering_scheme.csv


Here are what the first few lines of this data frame look like:


```python
display(HTML(aafreqs_df.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>year</th>
      <th>isite</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>natural_count</th>
      <th>natural_frequency</th>
      <th>site</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2017.0</td>
      <td>0</td>
      <td>M</td>
      <td>M</td>
      <td>731</td>
      <td>1.0</td>
      <td>-16</td>
    </tr>
    <tr>
      <td>2014.0</td>
      <td>0</td>
      <td>M</td>
      <td>M</td>
      <td>161</td>
      <td>1.0</td>
      <td>-16</td>
    </tr>
    <tr>
      <td>2017.5</td>
      <td>0</td>
      <td>M</td>
      <td>M</td>
      <td>187</td>
      <td>1.0</td>
      <td>-16</td>
    </tr>
    <tr>
      <td>2013.0</td>
      <td>0</td>
      <td>M</td>
      <td>M</td>
      <td>500</td>
      <td>1.0</td>
      <td>-16</td>
    </tr>
    <tr>
      <td>2011.0</td>
      <td>0</td>
      <td>M</td>
      <td>M</td>
      <td>260</td>
      <td>1.0</td>
      <td>-16</td>
    </tr>
  </tbody>
</table>


Now write this amino-acid frequency data frame to a file:


```python
aafreqs_file = os.path.join(resultsdir, 'aafreqs.csv')
print(f"Writing amino-acid frequencies to {aafreqs_file}")
aafreqs_df.to_csv(aafreqs_file, index=False)
```

    Writing amino-acid frequencies to results/natseqs/aafreqs.csv


## Examine sites of strong immune selection
Now we examine the amino-acid frequency changes at sites of strong immune selection.

First, we merge the amino-acid frequencies with the data frame containing the immune-selection values.
We keep only sites of major selection from contemporary human sera, which we define as:

 - being among the "contemporaneous" human sera, which is the "VIDD sera" set.
 - being "significant" sites according to the criterion used when analyzing the mutational antigenic profiling, which means that they have a value of `True` in the *zoom_site* column.
 
We also only keep amino acids at each site that have frequencies that exceed 5% in at least one year, and the re-normalize the frequencies for each site / year to sum to one:


```python
sel_df_file = os.path.join(config['avgdiffseldir'], 'avg_sel_tidy.csv')
print(f"Merging with selection values in {sel_df_file}")
sel_df = pd.read_csv(sel_df_file, low_memory=False)

# add immune selection values and get sites of major selection
df = (pd.merge(aafreqs_df, sel_df)
      .query('serum_group == "VIDD_sera"')
      .query('zoom_site')
      [['year', 'site', 'mutation', 'natural_frequency']]
      .drop_duplicates()
      .rename(columns={'mutation': 'amino acid'})
      )

# subset to amino acids with high natural frequencies in at least one year
df = (df.merge(df
               .groupby(['site', 'amino acid'])
               ['natural_frequency']
               .max()
               .rename('max_freq')
               .reset_index()
               )
      .query('max_freq > 0.05')
      .drop(columns='max_freq')
      )

# re-normalize frequencies
df = (df.merge(df
               .groupby(['site', 'year'])
               ['natural_frequency']
               .sum()
               .rename('tot_freq')
               .reset_index()
               )
      .assign(natural_frequency=lambda x: x['natural_frequency'] / x['tot_freq'])
      .drop(columns='tot_freq')
      )
```

    Merging with selection values in results/avgdiffsel/avg_sel_tidy.csv


Plot frequencies per year for these sites:


```python
p = (
 ggplot(df, aes('year', 'natural_frequency', color='amino acid')) +
 geom_line() +
 facet_wrap('site', nrow=2) +
 theme(figure_size=(9, 4))
 )

_ = p.draw()
plotfile = os.path.join(resultsdir, 'natural_freqs_plot.pdf')
p.save(plotfile)
```


![png](analyze_natseqs_files/analyze_natseqs_29_0.png)



```python

```
