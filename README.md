# Neff Calculation Script

This script calculates the effective number of sequences (Neff) from a multiple sequence alignment (MSA) file using the following equation:

$$ Neff = \sum_{i=1}^{N} {1 \over weight_{i}} $$

where $N$ is the number of sequences in the MSA. For each sequence $i$, $weight_{i}$ is the sum of $weight_{ij}$ over all sequences in the MSA, while $weight_{ij}$ is the sequence identity between any homologous sequence $i$ and $j$ in the MSA. The minimum sequence identity between $i$ and $j$ is set to 80% by default but can be defined by the user (see below).

## Requirements

- Python 3
- numpy

## Arguments

- `--msa_file`: Path to an MSA file. Supported formats are `.a3m`, `.fasta`, `.sto`, and `.hhr`.
- `--af_dir`: Path to an AlphaFold2 output directory. This will automatically find the MSA files and calculate the total number of sequences/Neff.
- `--cf_dir`: Path to a ColabFold output directory. This will automatically find the MSA file and calculate the total number of sequences/Neff.
- `--seq_id`: Optional: Sequence identity threshold for calculating Neff. The default is `0.8` in line with the 80% value used in the AlphaFold2 paper ([Jumper *et al.* 2021](https://doi.org/10.1038/s41586-021-03819-2)). However, be aware that other thresholds are also used: e.g. `0.62` ([Wu *et al.* 2020](https://doi.org/10.1093/bioinformatics/btz679)).
- `--input_seq`: Optional: The input sequence used to calculate the MSA (fasta format). This is used to calculate the Neff/length.

## Usage

```bash
python get_neff.py --msa_file <path_to_msa_file> [--seq_id <sequence_identity_threshold> --input_seq <path_to_input_sequence>]
python get_neff.py --af_dir <path_to_af_dir> [--seq_id <sequence_identity_threshold> --input_seq <path_to_input_sequence>]
python get_neff.py --cf_dir <path_to_cf_dir> [--seq_id <sequence_identity_threshold> --input_seq <path_to_input_sequence>]
```

Note: Can be a little slow on large MSAs (>10000 sequences)

## Output

Results are printed to the terminal in the format:
```bash
Neff = <neff value>
Neff/len = <neff/len value>
Total sequences = <total_sequences>
```
The results are also stored in an output file called `neff.csv`.

## Acknowledgments

This script uses lightly modified code from [AlphaFold](https://github.com/google-deepmind/alphafold) and [ColabDesign](https://github.com/sokrypton/ColabDesign/tree/ed4b01354928b60cd1347f570e9b248f78f11c6d).

