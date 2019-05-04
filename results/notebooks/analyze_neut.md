
<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Fit-and-plot-neutralization-curves" data-toc-modified-id="Fit-and-plot-neutralization-curves-1">Fit and plot neutralization curves</a></span><ul class="toc-item"><li><span><a href="#Import-Python-modules-/-packages" data-toc-modified-id="Import-Python-modules-/-packages-1.1">Import Python modules / packages</a></span></li><li><span><a href="#Configuration-and-setup" data-toc-modified-id="Configuration-and-setup-1.2">Configuration and setup</a></span></li><li><span><a href="#Read-neutralization-data" data-toc-modified-id="Read-neutralization-data-1.3">Read neutralization data</a></span></li><li><span><a href="#Fit-and-plot-all-neutralization-curves" data-toc-modified-id="Fit-and-plot-all-neutralization-curves-1.4">Fit and plot all neutralization curves</a></span></li><li><span><a href="#Colored-plots-for-figures" data-toc-modified-id="Colored-plots-for-figures-1.5">Colored plots for figures</a></span></li></ul></li></ul></div>

# Fit and plot neutralization curves
In this notebook we will plot neutralization curves from GFP-based neutralization assays. 
The GFP-based neutralization assay system is described in detail [here](https://github.com/jbloomlab/flu_PB1flank-GFP_neut_assay).

All the curves plotted here represent the mean and standard deviation of three replicates, with each replicate in a separate column of a 96-well plate.

The curves fit are the Hill-style neutralization functions fit and plotted by the [neutcurve](https://jbloomlab.github.io/neutcurve/) package.

## Import Python modules / packages


```python
import os
import warnings

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import pandas as pd

import yaml

import neutcurve
from neutcurve.colorschemes import CBPALETTE
import neutcurve.parse_excel

print(f"Using neutcurve version {neutcurve.__version__}")
```

    Using neutcurve version 0.2.dev0


Suppress warnings that can clutter output:


```python
warnings.simplefilter('ignore')
```

## Configuration and setup
Read general configuration from [config.yaml](config.yaml), which in turn specifies additional configuration files:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Read the neutralization assay configuration from the specified file:


```python
print(f"Reading neutralization assay setup from {config['neut_config']}")

with open(config['neut_config']) as f:
    neut_config = yaml.safe_load(f)
```

    Reading neutralization assay setup from data/neut_assays/neut_config.yaml


Get the output directory:


```python
outdir = config['neutresultsdir']
os.makedirs(outdir, exist_ok=True)
print(f"Output will be written to {outdir}")
```

    Output will be written to results/neutralization_assays


## Read neutralization data

Next, for each dict in *neut_config*, we use [neutcurve.parse_excel.parseRachelStyle2019](https://jbloomlab.github.io/neutcurve/neutcurve.parse_excel.html#neutcurve.parse_excel.parseRachelStyle2019) to parse the raw Excel files to create a tidy data frame appropriate for passing to [neutcurve.CurveFits](https://jbloomlab.github.io/neutcurve/neutcurve.curvefits.html#neutcurve.curvefits.CurveFits). 
We then concatenate all the tidy data frames to get our neutralization data.
This is essentially the workflow [explained here](https://jbloomlab.github.io/neutcurve/rachelstyle2019_example.html):


```python
neutdata = []  # store all data frame, then concatenate at end

for sampledict in neut_config:
    assert len(sampledict) == 1
    sampleset, kwargs = list(sampledict.items())[0]
    print(f"Parsing data for {sampleset}...")
    neutdata.append(neutcurve.parse_excel.parseRachelStyle2019(**kwargs))

neutdata = pd.concat(neutdata)
print(f"Read data for {len(neutdata.groupby('serum'))} sera and "
      f"{len(neutdata.groupby(['serum', 'virus']))} serum / virus pairs.")
```

    Parsing data for VIDD1...
    Parsing data for VIDD2...
    Parsing data for VIDD3...
    Parsing data for VIDD4...
    Parsing data for VIDD5...
    Parsing data for 557v1...
    Parsing data for 557v2...
    Parsing data for 574v1...
    Parsing data for 574v2...
    Parsing data for 589v1...
    Parsing data for 589v2...
    Parsing data for 571v1...
    Parsing data for 571v2...
    Parsing data for ferret-Pitt-1-preinf...
    Parsing data for ferret-Pitt-1-postinf...
    Parsing data for ferret-Pitt-2-preinf...
    Parsing data for ferret-Pitt-2-postinf...
    Parsing data for ferret-Pitt-3-preinf...
    Parsing data for ferret-Pitt-3-postinf...
    Parsing data for ferret-WHO...
    Parsing data for ferret-WHO-Victoria2011...
    Parsing data for antibody-5A01...
    Parsing data for antibody-3C06...
    Parsing data for antibody-3C04...
    Parsing data for antibody-4C01...
    Parsing data for antibody-4F03...
    Parsing data for antibody-1C04...
    Read data for 27 sera and 91 serum / virus pairs.


These data look like this:


```python
display(HTML(neutdata.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2010-age-21</td>
      <td>wt</td>
      <td>1</td>
      <td>0.000024</td>
      <td>1.010575</td>
    </tr>
    <tr>
      <td>2010-age-21</td>
      <td>wt</td>
      <td>1</td>
      <td>0.000046</td>
      <td>0.981598</td>
    </tr>
    <tr>
      <td>2010-age-21</td>
      <td>wt</td>
      <td>1</td>
      <td>0.000086</td>
      <td>0.997023</td>
    </tr>
    <tr>
      <td>2010-age-21</td>
      <td>wt</td>
      <td>1</td>
      <td>0.000163</td>
      <td>0.954407</td>
    </tr>
    <tr>
      <td>2010-age-21</td>
      <td>wt</td>
      <td>1</td>
      <td>0.000308</td>
      <td>1.005039</td>
    </tr>
  </tbody>
</table>


Write the neutralization data to a CSV file in our output directory:


```python
neutdatafile = os.path.join(outdir, 'neutdata.csv')
neutdata.to_csv(neutdatafile, index=False)
print(f"Wrote neutralization data to {neutdatafile}")
```

    Wrote neutralization data to results/neutralization_assays/neutdata.csv


## Fit and plot all neutralization curves

Now we fit the neutralization curves with a neutcurve.CurveFits:


```python
fits = neutcurve.CurveFits(neutdata)
```

Make a big panel of plots of the across-replicate averages for all sera and antibodies:


```python
for is_antibody, ptype, xlabel in [(False, 'sera', 'serum dilution'),
                                   (True, 'antibody', 'concentration ($\mu$g/ml)')]:
    plotfile = os.path.join(outdir, f"all_{ptype}_plots.pdf")
    print(f"\nPlotting all {ptype} curves to {plotfile}...")
    fig, _ = fits.plotSera(
                sera=[s for s in fits.sera if ('antibody' in s) == is_antibody],
                xlabel=xlabel,
                max_viruses_per_subplot=6,
                )
    display(fig)
    fig.savefig(plotfile)
    plt.close(fig)
```

    
    Plotting all sera curves to results/neutralization_assays/all_sera_plots.pdf...



![png](analyze_neut_files/analyze_neut_21_1.png)


    
    Plotting all antibody curves to results/neutralization_assays/all_antibody_plots.pdf...



![png](analyze_neut_files/analyze_neut_21_3.png)


Now get the curve fit parameters (e.g., IC50s):


```python
fitparams = fits.fitParams()
```

Here is a list of all of the IC50s for the replicate-averages for each serum / virus pair.
Note that for sera there are dilutions, and for antibodies that are $\mu$g/ml:


```python
display(HTML(fitparams
             .query('replicate == "average"')
             .drop(columns=['midpoint', 'top', 'bottom', 'replicate',
                            'ic50', 'ic50_bound', 'slope', 'nreplicates'])
             .to_html(index=False)
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>ic50_str</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2010-age-21</td>
      <td>wt</td>
      <td>0.00118</td>
    </tr>
    <tr>
      <td>2010-age-21</td>
      <td>F193D</td>
      <td>&gt;0.014</td>
    </tr>
    <tr>
      <td>2010-age-21</td>
      <td>syn</td>
      <td>0.00164</td>
    </tr>
    <tr>
      <td>2010-age-21</td>
      <td>K189D</td>
      <td>0.0018</td>
    </tr>
    <tr>
      <td>2010-age-21</td>
      <td>L157D</td>
      <td>0.00107</td>
    </tr>
    <tr>
      <td>2009-age-53a</td>
      <td>wt</td>
      <td>0.00124</td>
    </tr>
    <tr>
      <td>2009-age-53a</td>
      <td>L157D</td>
      <td>0.00481</td>
    </tr>
    <tr>
      <td>2009-age-53a</td>
      <td>K160T</td>
      <td>0.00338</td>
    </tr>
    <tr>
      <td>2009-age-53a</td>
      <td>F193D</td>
      <td>0.00234</td>
    </tr>
    <tr>
      <td>2009-age-53a</td>
      <td>syn</td>
      <td>0.00153</td>
    </tr>
    <tr>
      <td>2009-age-53b</td>
      <td>wt</td>
      <td>0.00172</td>
    </tr>
    <tr>
      <td>2009-age-53b</td>
      <td>L157D</td>
      <td>0.00931</td>
    </tr>
    <tr>
      <td>2009-age-53b</td>
      <td>K160T</td>
      <td>0.00387</td>
    </tr>
    <tr>
      <td>2009-age-53b</td>
      <td>F193D</td>
      <td>0.00418</td>
    </tr>
    <tr>
      <td>2009-age-53b</td>
      <td>F159G</td>
      <td>0.00174</td>
    </tr>
    <tr>
      <td>2009-age-64</td>
      <td>wt</td>
      <td>0.000318</td>
    </tr>
    <tr>
      <td>2009-age-64</td>
      <td>F159G</td>
      <td>0.0108</td>
    </tr>
    <tr>
      <td>2009-age-64</td>
      <td>K189D</td>
      <td>0.000258</td>
    </tr>
    <tr>
      <td>2009-age-64</td>
      <td>F193D</td>
      <td>0.000386</td>
    </tr>
    <tr>
      <td>2009-age-65</td>
      <td>wt</td>
      <td>0.000328</td>
    </tr>
    <tr>
      <td>2009-age-65</td>
      <td>F193D</td>
      <td>&gt;0.0109</td>
    </tr>
    <tr>
      <td>2009-age-65</td>
      <td>K160T</td>
      <td>0.00145</td>
    </tr>
    <tr>
      <td>2009-age-65</td>
      <td>syn</td>
      <td>0.000412</td>
    </tr>
    <tr>
      <td>2009-age-65</td>
      <td>N121E</td>
      <td>0.000332</td>
    </tr>
    <tr>
      <td>2009-age-65</td>
      <td>F159G</td>
      <td>0.000805</td>
    </tr>
    <tr>
      <td>2015-age-25-prevacc</td>
      <td>wt</td>
      <td>0.00126</td>
    </tr>
    <tr>
      <td>2015-age-25-prevacc</td>
      <td>F159G</td>
      <td>&gt;0.0102</td>
    </tr>
    <tr>
      <td>2015-age-25-prevacc</td>
      <td>R220D</td>
      <td>&gt;0.0102</td>
    </tr>
    <tr>
      <td>2015-age-25-prevacc</td>
      <td>K189D</td>
      <td>0.000868</td>
    </tr>
    <tr>
      <td>2015-age-25-vacc</td>
      <td>wt</td>
      <td>8.54e-05</td>
    </tr>
    <tr>
      <td>2015-age-25-vacc</td>
      <td>F159G</td>
      <td>&gt;0.00117</td>
    </tr>
    <tr>
      <td>2015-age-25-vacc</td>
      <td>syn</td>
      <td>0.000141</td>
    </tr>
    <tr>
      <td>2015-age-29-prevacc</td>
      <td>wt</td>
      <td>0.00637</td>
    </tr>
    <tr>
      <td>2015-age-29-prevacc</td>
      <td>F159G</td>
      <td>0.00954</td>
    </tr>
    <tr>
      <td>2015-age-29-prevacc</td>
      <td>K144E</td>
      <td>0.00561</td>
    </tr>
    <tr>
      <td>2015-age-29-vacc</td>
      <td>wt</td>
      <td>0.000299</td>
    </tr>
    <tr>
      <td>2015-age-29-vacc</td>
      <td>K144E</td>
      <td>0.000559</td>
    </tr>
    <tr>
      <td>2015-age-29-vacc</td>
      <td>F159G</td>
      <td>0.000293</td>
    </tr>
    <tr>
      <td>2015-age-48-prevacc</td>
      <td>wt</td>
      <td>&gt;0.00617</td>
    </tr>
    <tr>
      <td>2015-age-48-vacc</td>
      <td>wt</td>
      <td>8.33e-05</td>
    </tr>
    <tr>
      <td>2015-age-48-vacc</td>
      <td>K189D</td>
      <td>0.000413</td>
    </tr>
    <tr>
      <td>2015-age-48-vacc</td>
      <td>F159G</td>
      <td>0.000149</td>
    </tr>
    <tr>
      <td>2015-age-48-vacc</td>
      <td>R220D</td>
      <td>5.11e-06</td>
    </tr>
    <tr>
      <td>2015-age-48-vacc</td>
      <td>K144E</td>
      <td>7.23e-05</td>
    </tr>
    <tr>
      <td>2015-age-49-prevacc</td>
      <td>wt</td>
      <td>&gt;0.00617</td>
    </tr>
    <tr>
      <td>2015-age-49-vacc</td>
      <td>wt</td>
      <td>0.00343</td>
    </tr>
    <tr>
      <td>2015-age-49-vacc</td>
      <td>G75(HA2)H</td>
      <td>0.00277</td>
    </tr>
    <tr>
      <td>ferret-Pitt-1-preinf</td>
      <td>wt</td>
      <td>&gt;0.00617</td>
    </tr>
    <tr>
      <td>ferret-Pitt-1-postinf</td>
      <td>wt</td>
      <td>9.84e-05</td>
    </tr>
    <tr>
      <td>ferret-Pitt-1-postinf</td>
      <td>K189D</td>
      <td>0.000292</td>
    </tr>
    <tr>
      <td>ferret-Pitt-1-postinf</td>
      <td>F193D</td>
      <td>0.000484</td>
    </tr>
    <tr>
      <td>ferret-Pitt-1-postinf</td>
      <td>syn</td>
      <td>9.9e-05</td>
    </tr>
    <tr>
      <td>ferret-Pitt-1-postinf</td>
      <td>F159G</td>
      <td>0.000186</td>
    </tr>
    <tr>
      <td>ferret-Pitt-1-postinf</td>
      <td>K144E</td>
      <td>0.000176</td>
    </tr>
    <tr>
      <td>ferret-Pitt-2-preinf</td>
      <td>wt</td>
      <td>&gt;0.00206</td>
    </tr>
    <tr>
      <td>ferret-Pitt-2-postinf</td>
      <td>wt</td>
      <td>0.000343</td>
    </tr>
    <tr>
      <td>ferret-Pitt-2-postinf</td>
      <td>K189D</td>
      <td>0.000724</td>
    </tr>
    <tr>
      <td>ferret-Pitt-2-postinf</td>
      <td>F193D</td>
      <td>0.0012</td>
    </tr>
    <tr>
      <td>ferret-Pitt-2-postinf</td>
      <td>K144E</td>
      <td>0.000636</td>
    </tr>
    <tr>
      <td>ferret-Pitt-3-preinf</td>
      <td>wt</td>
      <td>&gt;0.00617</td>
    </tr>
    <tr>
      <td>ferret-Pitt-3-postinf</td>
      <td>wt</td>
      <td>0.000375</td>
    </tr>
    <tr>
      <td>ferret-Pitt-3-postinf</td>
      <td>K189D</td>
      <td>0.00113</td>
    </tr>
    <tr>
      <td>ferret-Pitt-3-postinf</td>
      <td>F193D</td>
      <td>0.00505</td>
    </tr>
    <tr>
      <td>ferret-Pitt-3-postinf</td>
      <td>K160T</td>
      <td>0.000471</td>
    </tr>
    <tr>
      <td>ferret-Pitt-3-postinf</td>
      <td>L157D</td>
      <td>0.000265</td>
    </tr>
    <tr>
      <td>ferret-WHO</td>
      <td>wt</td>
      <td>0.00476</td>
    </tr>
    <tr>
      <td>ferret-WHO-Victoria2011</td>
      <td>wt</td>
      <td>0.00192</td>
    </tr>
    <tr>
      <td>ferret-WHO-Victoria2011</td>
      <td>K189D</td>
      <td>0.0115</td>
    </tr>
    <tr>
      <td>ferret-WHO-Victoria2011</td>
      <td>F193D</td>
      <td>&gt;0.0125</td>
    </tr>
    <tr>
      <td>ferret-WHO-Victoria2011</td>
      <td>F159G</td>
      <td>0.00601</td>
    </tr>
    <tr>
      <td>ferret-WHO-Victoria2011</td>
      <td>K160T</td>
      <td>0.00871</td>
    </tr>
    <tr>
      <td>antibody-5A01</td>
      <td>wt</td>
      <td>0.0884</td>
    </tr>
    <tr>
      <td>antibody-5A01</td>
      <td>F159G</td>
      <td>&gt;1</td>
    </tr>
    <tr>
      <td>antibody-5A01</td>
      <td>syn</td>
      <td>0.142</td>
    </tr>
    <tr>
      <td>antibody-5A01</td>
      <td>K160T</td>
      <td>0.549</td>
    </tr>
    <tr>
      <td>antibody-3C06</td>
      <td>wt</td>
      <td>0.00726</td>
    </tr>
    <tr>
      <td>antibody-3C06</td>
      <td>F159G</td>
      <td>&gt;0.5</td>
    </tr>
    <tr>
      <td>antibody-3C06</td>
      <td>K160T</td>
      <td>0.0269</td>
    </tr>
    <tr>
      <td>antibody-3C04</td>
      <td>wt</td>
      <td>0.00774</td>
    </tr>
    <tr>
      <td>antibody-3C04</td>
      <td>K160T</td>
      <td>&gt;0.5</td>
    </tr>
    <tr>
      <td>antibody-3C04</td>
      <td>K160S</td>
      <td>0.333</td>
    </tr>
    <tr>
      <td>antibody-4C01</td>
      <td>wt</td>
      <td>0.0207</td>
    </tr>
    <tr>
      <td>antibody-4C01</td>
      <td>F193D</td>
      <td>&gt;0.5</td>
    </tr>
    <tr>
      <td>antibody-4C01</td>
      <td>syn</td>
      <td>0.0112</td>
    </tr>
    <tr>
      <td>antibody-4C01</td>
      <td>K160T</td>
      <td>0.0202</td>
    </tr>
    <tr>
      <td>antibody-4F03</td>
      <td>wt</td>
      <td>0.362</td>
    </tr>
    <tr>
      <td>antibody-4F03</td>
      <td>N121E</td>
      <td>&gt;3</td>
    </tr>
    <tr>
      <td>antibody-4F03</td>
      <td>F193D</td>
      <td>0.13</td>
    </tr>
    <tr>
      <td>antibody-4F03</td>
      <td>F159G</td>
      <td>0.0649</td>
    </tr>
    <tr>
      <td>antibody-1C04</td>
      <td>wt</td>
      <td>0.0868</td>
    </tr>
    <tr>
      <td>antibody-1C04</td>
      <td>K82A</td>
      <td>&gt;5</td>
    </tr>
  </tbody>
</table>


Write all of the fit parameters (including IC50s) for all replicates to a file:


```python
fitfile = os.path.join(outdir, 'fitparams.csv')
print(f"Writing fit parameters to {fitfile}")
fitparams.to_csv(fitfile, index=False, float_format='%.3g')
```

    Writing fit parameters to results/neutralization_assays/fitparams.csv


## Colored plots for figures
Get the colors that we have specified for viruses in each serum group:


```python
with open(config['mutation_colors_and_markers']) as f:
    virus_to_color_marker = yaml.safe_load(f)
```

Specify the extensions for the figure files:


```python
fig_extensions = ['.svg', '.pdf']
```

Write figures to this directory:


```python
figsdir = config['figsdir']
os.makedirs(figsdir, exist_ok=True)
```

Get data frame with all sera:


```python
with open(config['serum_info']) as f:
    sera_df = (pd.DataFrame.from_dict(yaml.safe_load(f))
               .transpose()
               .rename_axis('serum')
               .reset_index()
               .assign(name=lambda x:
                            x['species'].map(lambda s: '' if pd.isnull(s) or
                                             s == 'human' else s + '-') +
                            x['name']
                       )
               )
```

Now draw the plots as single-column figures so they can be combined with the logo plots:


```python
for serum_group, df in sera_df.groupby('group'):
    if serum_group not in virus_to_color_marker:
        continue
    sera = df['name'].unique()
    if 'antibody' in serum_group:
        xlabel = 'concentration ($\mu$g/ml)'
    else:
        xlabel = 'serum dilution'
    colors = virus_to_color_marker[serum_group]['colors']
    print(f"\n**************** {serum_group} ****************")
    fig, _ = fits.plotSera(
                sera=df['name'].unique(),
                viruses=colors.keys(),
                virus_to_color_marker=colors,
                xlabel=xlabel,
                max_viruses_per_subplot=len(colors),
                ncol=1,
                )
    display(fig)
    for ext in fig_extensions:
        plotfile = os.path.join(config['figsdir'], f"{serum_group}_neut{ext}")
        print(f"Saving to {plotfile}")
        fig.savefig(plotfile)
    plt.close(fig)
```

    
    **************** Hensley_sera ****************



![png](analyze_neut_files/analyze_neut_37_1.png)


    Saving to results/figures/Hensley_sera_neut.svg
    Saving to results/figures/Hensley_sera_neut.pdf
    
    **************** VIDD_sera ****************



![png](analyze_neut_files/analyze_neut_37_3.png)


    Saving to results/figures/VIDD_sera_neut.svg
    Saving to results/figures/VIDD_sera_neut.pdf
    
    **************** antibody_lower_head ****************



![png](analyze_neut_files/analyze_neut_37_5.png)


    Saving to results/figures/antibody_lower_head_neut.svg
    Saving to results/figures/antibody_lower_head_neut.pdf
    
    **************** antibody_region_B ****************



![png](analyze_neut_files/analyze_neut_37_7.png)


    Saving to results/figures/antibody_region_B_neut.svg
    Saving to results/figures/antibody_region_B_neut.pdf
    
    **************** ferret ****************



![png](analyze_neut_files/analyze_neut_37_9.png)


    Saving to results/figures/ferret_neut.svg
    Saving to results/figures/ferret_neut.pdf



```python

```
