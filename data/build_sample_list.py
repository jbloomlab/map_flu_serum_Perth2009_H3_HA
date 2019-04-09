"""Build list of samples from Juhye's experiment list."""


import os
import glob

import numpy
import pandas as pd


def main():

    # map Juhye's FASTQ keys to directories
    paths = {'A':'/shared/ngs/illumina/jmlee34/'
             '190326_D00300_0706_AHWYTCBCX2/Unaligned/'
             'Project_jmlee34/',
             'B':'/shared/ngs/illumina/jmlee34/'
                 '190116_D00300_0668_AHVHCFBCX2/Unaligned/'
                 'Project_jmlee34/'
             }

    # map Juhye's FASTQ keys to sequencing dates
    dates = {'A':'2019-03-26',
             'B':'2019-01-16',
             }

    # map Juhye's library notations to library names
    # note we use libraries starting 1, 2, ... rather than
    # the experiment numbering as this makes more sense
    # to outsider if libraries start at 1
    libraries = {'L4':'lib1',
                 'Lib4':'lib1',
                 'L5':'lib2',
                 'Lib5':'lib2', 
                 'L6':'lib3', 
                 'Lib6':'lib3'
                 }

    # read in Juhye's experiment list, add R1, date, library
    exp_list_file = './data/experiment_list.csv'
    print(f"Reading experiment list from {exp_list_file}")
    exp_list = (
        pd.read_csv(exp_list_file)
        .assign(
            fastq_path=lambda x: x.FASTQkey.map(paths),
            R1=lambda x: x.fastq_path + 'Sample_' + x.experiment +
                             '/*R1*.fastq.gz',
            n_R1files=lambda x: x.R1.apply(glob.glob).apply(len),
            date=lambda x: x.FASTQkey.map(dates),
            library=lambda x: x['name'].str
                              .extract('^(L(?:ib){0,1}\d)', expand=False)
                              .map(libraries),
            serum=lambda x: numpy.where(
                              x['name'].str.contains('plasmid'), 'plasmid',
                              numpy.where(
                                x['name'].str.contains('mock'), 'mock',
                                x['antibody']
                                )
                              ),
            serum_dilution=lambda x: x['serum_dilution']
                           .apply(round, ndigits=7)
                           .where((x['serum'] != 'plasmid') &
                                  (x['serum'] != 'mock')),
            percent_infectivity=lambda x: x['percent_infectivity'].where(
                                x['serum'] != 'plasmid') 
            )
        .rename(columns={'name':'sample'})
        )

    # make sure we find R1 files for all samples
    if any(exp_list.n_R1files < 1):
        raise ValueError('Cannot find the following R1 files:\n' +
                         exp_list
                         .query('n_R1files < 1')
                         [['name', 'R1']]
                         .to_csv(index=False)
                         )

    # check library assigned for all samples not plasmid
    no_lib = exp_list.query('library.isnull() & serum != "plasmid"')
    if len(no_lib):
        raise ValueError('Failed to assign library for these samples:\n' +
                         no_lib[['name', 'library']].to_csv(index=False))

    # write sample list
    sample_list_file = 'sample_list.csv'
    print(f"Writing sample list to {sample_list_file}")
    (exp_list
     [['sample', 'serum', 'library', 'date',
       'serum_dilution', 'percent_infectivity', 'R1']]
     .to_csv('sample_list.csv', index=False)
     )


if __name__ == '__main__':
    main()
