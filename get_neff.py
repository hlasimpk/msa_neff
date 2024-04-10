import collections
import csv
import glob
import numpy as np
import os
import re
import string


ALPHABET = list("ARNDCQEGHILKMFPSTWYV-")

def parse_stockholm(filename):
    """Parse a Stockholm file and return a list of sequences."""

    with open(filename, 'r') as f_in:
        name_to_sequence = collections.OrderedDict()
        for line in f_in:
            line = line.strip()
            if not line or line.startswith(('#', '//')):
                continue
            name, sequence = line.split()
            if name not in name_to_sequence:
                name_to_sequence[name] = ''
            name_to_sequence[name] += sequence
            
    sequences = []
    keep_columns = []
    for seq_index, sequence in enumerate(name_to_sequence.values()):
        if seq_index == 0:
            query = sequence
            keep_columns = [i for i, res in enumerate(query) if res != '-']

        aligned_sequence = ''.join([sequence[c] for c in keep_columns])

        sequences.append(aligned_sequence)

    return sequences

def parse_fasta_a3m(filename, a3m=True, stop=100000):
  """Parse a fasta or a3m file and return a list of sequences."""

  if a3m:
    # for a3m files the lowercase letters are removed
    # as these do not align to the query sequence
    rm_lc = str.maketrans(dict.fromkeys(string.ascii_lowercase))

  header, sequence = [],[]
  lines = open(filename, "r")
  for line in lines:
    line = line.rstrip()
    if len(line) > 0:
      if line[0] == '#':
          continue
      if line[0] == ">":
        if len(header) == stop:
          break
        else:
          header.append(line[1:])
          sequence.append([])
      else:
        if a3m: 
            line = line.translate(rm_lc)
        else: 
            line = line.upper()
        sequence[-1].append(line)
  lines.close()
  sequence = [''.join(seq) for seq in sequence]

  return sequence

def parse_hhr(filename):
    """Parse a HHsearch HHR file and return a list of sequences."""
    with open(filename, 'r') as f_in:
        data = f_in.read()

    lines = data.splitlines()
    for line in lines:
        if line.startswith('Match_columns'):
            global seq_length 
            seq_length = int(line.split()[1].strip())

    block_starts = [i for i, line in enumerate(lines) if line.startswith('No ')]

    sequences = []
    if block_starts:
        block_starts.append(len(lines))  # Add the end of the final block.
        for i in range(len(block_starts) - 1):
            sequence = ['-' for i in range(seq_length)]
            for line in lines[block_starts[i]:block_starts[i + 1]][3:]:
                if line.startswith('Q ') and not line.startswith('Q Consensus'):
                    gaps = line.split()[3].count('-')
                    start = (int(line.split()[2]) - 1) - gaps
                    end = (int(line.split()[4]))

                if line.startswith('T '):
                    if (not line.startswith('T ss_dssp') and
                        not line.startswith('T ss_pred') and
                        not line.startswith('T Consensus')):
                        patt = r'[\t ]*([0-9]*) ([A-Z-]*)[\t ]*[0-9]* \([0-9]*\)'
                        groups = re.match(patt, line[17:]).groups()
                        delta_hit_sequence = groups[1]
                        for idx, res in enumerate(delta_hit_sequence):
                            sequence[start + idx] = res
            sequences.append(''.join(sequence))
    return sequences

def mk_msa(seqs):
  """one hot encode msa"""
  states = len(ALPHABET)
  a2n = {a:n for n,a in enumerate(ALPHABET)}
  msa_ori = np.array([[a2n.get(aa, states-1) for aa in seq] for seq in seqs])
  return np.eye(states)[msa_ori]


def calculate_neff(seq, seq_id=0.8):
    """Calculate the effective number of sequences."""
    msa = mk_msa(seq)
    msa_ident = np.tensordot(msa, msa, [[1,2],[1,2]]) / msa.shape[1]
    eff = 1 / (msa_ident >= seq_id).sum(-1)
    neff = eff.sum()
    return neff


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Calculate Neff')
    parser.add_argument('--input_seq', type=str, help='Fasta containing input sequence (used to calculate Neff/len)')
    parser.add_argument('--af_dir', type=str, help='AlphaFold output directory')
    parser.add_argument('--cf_dir', type=str, help='ColabFold output directory')
    parser.add_argument('--msa_file', type=str, help='MSA file (supported formates: a3m, fasta, sto, hhr)')
    parser.add_argument('--seq_id', type=float, default=0.62, help='Sequence identity threshold')
    args = parser.parse_args()

    if args.msa_file:
        _, ext = os.path.splitext(args.msa_file)
        if ext == '.a3m':
            sequences = parse_fasta_a3m(args.msa_file)
        elif ext == '.fasta':
            sequences = parse_fasta_a3m(args.msa_file, a3m=False)
        elif ext == '.sto':
            sequences = parse_stockholm(args.msa_file)
        elif ext == '.hhr':
            sequences = parse_hhr(args.msa_file)
        else:
            raise ValueError('Unsupported MSA format')
    elif args.af_dir:
        af2_msas = glob.glob(os.path.join(args.af_dir, 'msas', '*'))
        sequences = []
        for msa in af2_msas:
            if 'a3m' in msa:
                sequences += parse_fasta_a3m(msa)
            if 'sto' in msa:
                sequences += parse_stockholm(msa)
            if 'hhr' in msa:
                sequences += parse_hhr(msa)
    elif args.cf_dir:
        cf_msas = glob.glob(os.path.join(args.cf_dir, '*.a3m'))
        sequences = parse_fasta_a3m(cf_msas[0])
        
    neff = calculate_neff(sequences, seq_id=args.seq_id)

    print(f"Neff = {np.round(neff, 2)}")
    if args.input_seq:
        with open(args.input_seq, 'r') as f:
            next(f)
            seq_length = len(f.read().replace('\n', ''))
        print(f"Neff/len = {np.round(neff/seq_length, 2)}")
    print(f"Total sequences = {len(sequences)}")

    with open('neff.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Neff", "Neff/len", "Total sequences"])
        writer.writerow([neff, neff/seq_length if args.input_seq else '', len(sequences)])
