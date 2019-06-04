"""Papermill parameterize notebooks to map selection on structure."""


import os

import pandas as pd

import papermill

import yaml


def main():
    """Main body of script."""

    with open('config.yaml') as f:
        config = yaml.safe_load(f)

    os.makedirs(config['notebookdir'], exist_ok=True)

    avg_sel_df = os.path.join(config['avgdiffseldir'], 'avg_sel_tidy.csv')
    if not os.path.isfile(avg_sel_df):
        raise ValueError(f"Cannot find `avg_sel_df` {avg_sel_df}")
    
    # modify avg_sel_df to add PDB chain / site numbers
    (pd.read_csv(avg_sel_df, low_memory=False)
     .drop(columns=['pdb_chain', 'pdb_site'], errors='ignore')
     .merge(pd.read_csv(config['site_to_pdb']),
            on='site', how='left', validate='many_to_one')
     .to_csv(avg_sel_df, index=False, float_format='%.5g')
     )

    with open(config['figure_config']) as f:
        figure_config = yaml.safe_load(f)

    for fig, figconfig in figure_config['figures'].items():
        if 'sera_for_struct' not in figconfig:
            continue
        sera = figconfig['sera_for_struct']
        if len(sera) == 0:
            continue
        struct_nb = os.path.join(config['notebookdir'],
                                 f"map_on_struct_{fig}.ipynb")
        data_relpath = os.path.relpath(avg_sel_df, config['notebookdir'])
        outdir_relpath = os.path.relpath(config['structsdir'],
                                         config['notebookdir'])
        panelfig_relpath = os.path.relpath(os.path.join(config['figsdir'],
                                                        f"{fig}_struct.png"),
                                           config['notebookdir'])
        papermill.execute.execute_notebook(
                input_path=config['map_on_struct_template'],
                output_path=struct_nb,
                parameters={'data_csv': data_relpath,
                            'query_str': f"serum_name_formatted in {sera}",
                            'facet_col': 'serum_name_formatted',
                            'pdb': config['pdb_id'],
                            'orientation': figure_config['struct_orientation'],
                            'outdir': outdir_relpath,
                            'panel_fig': panelfig_relpath,
                            'panel_height': figure_config['struct_panel_height'],
                            },
                prepare_only=True,
                cwd=config['notebookdir'],
                )
        print(f"Created {struct_nb}")


if __name__ == '__main__':
    main()
