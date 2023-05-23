#!/usr/bin/env python3

"""
Created:      15/04/2023
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2023 C.A. Warmerdam

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import os
import sys
import argparse
import gzip

import numpy as np
import pandas as pd

# Metadata
__program__ = "eQTL sample overlap"
__author__ = "C.A. (Robert) Warmerdam"
__email__ = "c.a.warmerdam@umcg.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


# Constants

# Classes


# Classes
class SampleOverlapCalculator:
    def __init__(self, inclusion_path, maf_table,
                 table_name="filter_logs_full.log",
                 variant_inclusion_format="%s_SnpsToInclude.txt",
                 gene_inclusion_format="%s_GenesToInclude.txt"):

        self.overview_df = pd.read_table(os.path.join(inclusion_path, table_name), index_col='Dataset')

        self.overview_df['gene_inclusion_path'] = (
            self.overview_df.index.map(lambda name: os.path.join(gene_inclusion_format % name)))

        self.gene_inclusion_df = self.load_inclusion_df('gene_inclusion_path')

    def load_inclusion_df(self, column):
        # Create an empty dictionary to store dataframes
        dfs = {}

        # Generate dict of inclusion paths
        inclusion_paths = self.overview_df[column].to_dict()

        # Iterate over files in the directory
        for cohort, filepath in inclusion_paths.items():

            df = pd.read_csv(filepath)

            # Set the index to the 'ID' column
            df.set_index('ID', inplace=True)

            # Add the dataframe to the dictionary
            dfs[cohort] = df

        # Merge the dataframes on their index
        merged_df = pd.concat(dfs.values(), axis=1, keys=dfs.keys())

        # Replace NaN values with 0 and convert to boolean
        merged_df = merged_df.fillna(0).astype(bool)

        return merged_df

    def calculate_sample_overlap(self):

        # For each gene calculate the total number of samples that are included.
        sample_size_per_gene = (
            (self.gene_inclusion_df * self.overview_df['N'])
            .sum(axis=1))

        # For each gene-gene combination, calculate the number of overlapping samples.
        single_cohort_overlap = np.einsum('ic,jc->cij', self.gene_inclusion_df.values, self.gene_inclusion_df.values)
        overlapping_samples = (
            (single_cohort_overlap * self.overview_df['N'].values[:, np.newaxis, np.newaxis])
            .sum(axis=0))

        # For each gene-gene combination, ...
        return overlapping_samples / sample_size_per_gene[:,np.newaxis]



# Functions

# Main


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    parser = argparse.ArgumentParser(description='Annotate significant eQTL results with variant and gene information')



    args = parser.parse_args(argv[1:])
    print(args)

    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
