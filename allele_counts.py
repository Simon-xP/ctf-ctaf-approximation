"""
Simon Hampton
June 11th, 2024

allele_counts.py
    This program performs mutation calling on combined BAM files, requiring 
    a reference genome and list of mutations present in the sample.
"""

import os
import pysam
import pandas as pd

bam_dir = "merged_bam"
output_dir = "merged_tsv"
os.makedirs(output_dir, exist_ok=True)

# List of BAM files
bam_files = [
    os.path.join(bam_dir, f) for f in os.listdir(bam_dir)
    if f.startswith("merged") and f.endswith(".bam")
]

reference_genome = "/home/simonh/bam/reference.fa"
mutations = "FB22AN0514-PAIR.ruo.tsv"
mutations_df = pd.read_csv(mutations, sep='\t')

for bam_file in bam_files:
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
   
    allele_counts = []

    # Filter for mutant mutations only
    mutant_mutations_df = mutations_df[mutations_df['Type'] == 'Mutant']

    # Iterate over each mutation in the mutant list
    for _, mutation in mutant_mutations_df.iterrows():
        split_values = mutation['Chr.start'].split(':')
        if len(split_values) == 2:
            chrom = split_values[0]  # extract chromosome
            pos = int(split_values[1])  # extract position, convert to int
            ref = mutation['Ref']
            alt = mutation['Alt']

            # Fetch pileup columns at the specified position
            for pileupcolumn in bam.pileup(chrom, pos-1, pos, truncate=True):
                if pileupcolumn.pos == pos-1:
                    ref_count = 0
                    alt_count = 0
                    depth = pileupcolumn.n

                    # Count reference and alternate allele occurrences
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            if base == ref:
                                ref_count += 1
                            elif base == alt:
                                alt_count += 1

                    # Store information in data list
                    allele_counts.append({
                        'File': os.path.basename(bam_file),
                        'Chromosome': chrom,
                        'Position': pos,
                        'Depth': depth,
                        'ref_count': ref_count,
                        'alt_count': alt_count,
                        'ref_seq': ref,
                        'alt_seq': alt,
                        'Info': mutation['AAChange']  # Example: AAChange field from your mutation data
                    })
        else:
            # does not match expected format
            print(f"Unexpected Format: {mutation['Chr.start']}")

    bam.close()

    # Save data list as TSV
    df = pd.DataFrame(allele_counts)
    tsv_filename = os.path.join(output_dir, f"mutations_analysis_{os.path.basename(bam_file).split('.')[0]}.tsv")
    df.to_csv(tsv_filename, sep='\t', index=False)
    print(f"Saved results to {tsv_filename}")

