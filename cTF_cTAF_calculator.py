"""
Simon Hampton
July 2nd, 2024

cTF_cTAF_calculator.py
    This program uses the Bayesian Distribution to approximate cTF and cTAF
    using a logarithmic scale.
    Requires proper mutation_list.xlsx file and setup of merged_tsv folder

"""

import os
import pandas as pd
import numpy as np
from scipy.stats import poisson, uniform


# computes the logarithmic likelihood function based on poisson distribution
def log_likelihood(cfDNA_counts, local_depths, tumor_alle_freqs, cTF):
    log_likelihood_value = 0.0
    for i in range(len(cfDNA_counts)):
        expected_count = local_depths[i] * cTF * tumor_alle_freqs[i]
        if expected_count > 0:
            log_likelihood_value += poisson.logpmf(cfDNA_counts[i], expected_count)
        else:
            log_likelihood_value += float('-inf')  # If expected count is zero, log likelihood is -inf
    return log_likelihood_value

# uniform distribution reflects a lack of prior knowledge
def log_prior(cTF):
    return uniform.logpdf(cTF, loc=0, scale=1)  # uniform prior between 0 and 1

# log posterior computed by adding likelihood and prior
def log_posterior(cfDNA_counts, local_depths, tumor_alle_freqs, cTF):
    return log_likelihood(cfDNA_counts, local_depths, tumor_alle_freqs, cTF) + log_prior(cTF)


def estimate_cTF_and_cTAF(cfDNA_counts, local_depths, tumor_alle_freqs):
    if not (len(cfDNA_counts) == len(local_depths) == len(tumor_alle_freqs)):
        raise ValueError("All input lists must have the same length.")

    cTF_values = np.linspace(0, 1, 1000)  # uses values between 0 and 1 as example cTF values
    log_posterior_values = np.array([log_posterior(cfDNA_counts, local_depths, tumor_alle_freqs, cTF) for cTF in cTF_values])

    max_log_posterior = np.max(log_posterior_values)
    if np.isneginf(max_log_posterior):
        raise ValueError("Log posterior values are all negative infinity. Check input data and calculations.")
    
    # converting from log to regular space, normalize
    posterior_values = np.exp(log_posterior_values - max_log_posterior)
    normalized_posterior_values = posterior_values / np.sum(posterior_values)

    # compute cumulative sum, find median
    cTF_cumulative = np.cumsum(normalized_posterior_values)
    cTF_estimate = np.interp(0.5, cTF_cumulative, cTF_values)
    
    # compute cTAF estimate
    median_tumor_alle_freq = np.median(tumor_alle_freqs)
    cTAF_estimate = cTF_estimate * median_tumor_alle_freq

    return cTF_estimate, cTAF_estimate



if __name__ == "__main__":
    tissue_df = pd.read_excel("mutation_list.xlsx")
    tissue_df = tissue_df[tissue_df['Type'] == 'Mutant']

    # Split the 'chr.start' column into 'Chromosome' and 'Position'
    tissue_df[['Chromosome', 'Position']] = tissue_df['Chr.start'].str.split(':', expand=True)
    tissue_df['Position'] = tissue_df['Position'].astype(int)

    results_data = pd.DataFrame(columns=['Expected cTF', 'Trial #', 'Calculated cTF', 'Calculated cTAF', 'Difference'])

    for percent in [1, 2, 5, 7, 10, 20, 30, 40, 50]:
        for trial in [1, 2, 3, 4, 5]:
            plasma_df = pd.read_csv(f"merged_tsv/mutations_analysis_merged_{percent}_{trial}.tsv", sep='\t')

            # Ensure columns are correctly interpreted as numeric
            plasma_df['Position'] = pd.to_numeric(plasma_df['Position'], errors='coerce')
            plasma_df['alt_count'] = pd.to_numeric(plasma_df['alt_count'], errors='coerce')
            plasma_df['Depth'] = pd.to_numeric(plasma_df['Depth'], errors='coerce')

            # Merge dataframes based on the 'Position' column
            merged_df = pd.merge(plasma_df, tissue_df, on='Position')

            cfDNA_counts = merged_df['alt_count'].values
            local_depths = merged_df['Depth'].values
            tumor_alle_freqs = merged_df['AF'].values

            cTF_estimate, cTAF_estimate = estimate_cTF_and_cTAF(cfDNA_counts, local_depths, tumor_alle_freqs)

            print(f"percent: {percent}")
            print(f"trial: {trial}")
            print(f"cTF: {cTF_estimate}")
            print(f"cTAF: {cTAF_estimate}")

            decimal = percent / 100
            difference = abs(cTF_estimate - decimal)

            # Create a new row as a DataFrame and append to results_data
            new_row = pd.DataFrame([[decimal, trial, cTF_estimate, cTAF_estimate, difference]],
                                   columns=['Expected cTF', 'Trial #', 'Calculated cTF', 'Calculated cTAF', 'Difference'])
            results_data = pd.concat([results_data, new_row], ignore_index=True)

    results_data.to_csv("cTF_cTAF_results.csv", index=False)
    print(f"Saved results to cTF_cTAF_results.csv")
