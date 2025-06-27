#!/usr/bin/env python

"""Usage:
    score.py [--fmin <fmin>] 
             [--fmax <fmax>]
             [--input <input_dir>]
             [--output <output_dir>]
             <sample_id>

Arguments:
    <sample_id>   The ID of the sample to process.

Options:
    --fmin <fmin>           The minimum fragment length [default: 100].
    --fmax <fmax>           The maximum fragment length [default: 500].
    --input <input_dir>     Path to input BAM files [default: /media/stagiaireNFS4/maya/INTENSO/PADA-1]
    --output <output_dir>   Path to save results [default: /media/stagiaireNFS4/maya/results]
"""

import numpy as np
import pandas as pd
import pysam
import os
from docopt import docopt

# Constants
bedpath = '/media/jnoirelNFS3/INTENSO/YODA_hg38.bed'
epsilon = 1e-6
hw = 60
window = 2 * hw + 1

def parse_chr(chr_num_list):
    chr_list = []
    for chr_num in chr_num_list.split(','):
        if '-' in chr_num:
            start, end = map(int, chr_num.split('-'))
            chr_list.extend([f'chr{x}' for x in range(start, end + 1)])
        else:
            chr_list.append('chr' + chr_num)
    return chr_list

def process_chr(bamfile, bed, sample, chrom, size, fmin, fmax, window):
    print(f'Processing {chrom}')
    print('...Initialisation')

    counts = {'sample': sample, 'chr': chrom,
              'allfrag': 0, 'allfrag-crit': 0,
              'nodup': 0, 'nodup-crit': 0}

    FD_a = np.zeros(size)
    FD_n = np.zeros(size)

    FS_a = np.zeros(size)
    FS_n = np.zeros(size)

    PS_plus = np.zeros(size)
    PS_minus = np.zeros(size)
    PS_plus_nodup = np.zeros(size)
    PS_minus_nodup = np.zeros(size)

    print('...Reading the reads')

    for read in bamfile.fetch(chrom):
        if read.is_reverse or not read.is_proper_pair:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        flen = read.template_length
        l = read.reference_start
        r = l + flen

        if r >= size or l < 0:
            continue

        counts['allfrag'] += 1
        if read.is_duplicate:
            counts['nodup'] += 1

        if not (fmin <= flen <= fmax):
            continue

        FD_a[l:r] += 1
        FS_a[l] += 1
        FS_a[r] += 1

        if flen < window:
            PS_minus[max(0, l - hw):min(size, r + hw)] += 1
        else:
            PS_plus[l + hw:r - hw] += 1
            PS_minus[max(0, l - hw):(l + hw)] += 1
            PS_minus[(r - hw):min(size, r + hw)] += 1

        counts['allfrag-crit'] += 1
        if read.is_duplicate:
            counts['nodup-crit'] += 1
            
        if read.is_duplicate:
           continue

        FD_n[l:r] += 1
        FS_n[l] += 1
        FS_n[r] += 1
        
        if flen < window:
           PS_minus_nodup[max(0, l - hw):min(size, r + hw)] += 1
        else:
           PS_plus_nodup[l + hw:r - hw] += 1
           PS_minus_nodup[max(0, l - hw):(l + hw)] += 1
           PS_minus_nodup[(r - hw):min(size, r + hw)] += 1

    print('...Wrapping up and intersecting with the BED file')

    dat = pd.DataFrame({
        'sample': sample,
        'chromosome': chrom,
        'position': np.arange(size, dtype=int),
        'frag_depth_all': FD_a,
        'frag_suscp_all': FS_a,
        'frag_depth_nodup': FD_n,
        'frag_suscp_nodup': FS_n,
        'protection_score_plus': PS_plus,
        'protection_score_minus': PS_minus,
        'protection_score_plus_nodup': PS_plus_nodup,
        'protection_score_minus_nodup': PS_minus_nodup,
    })

    bed_chrom = bed[bed['chr'] == chrom]
    results = []

    for _, row in bed_chrom.iterrows():
        region_data = dat[(dat['position'] >= row['start']) & (dat['position'] <= row['end'])]
        if not region_data.empty:
            region_data = region_data.assign(region_id=row['region_id'])
            results.append(region_data)

    if results:
        return pd.concat(results, ignore_index=True)
    else:
        return pd.DataFrame(columns=dat.columns)

def main(sample, chrs, bed, fmin, fmax, window, input_dir, output_dir):
    bampath = os.path.join(input_dir, sample, f'{sample}_dup.bam')
    results = []
    print('SCORE: Reading the BAM file')

    with pysam.AlignmentFile(bampath, 'rb') as bamfile:
        sizes = dict(zip(bamfile.references, bamfile.lengths))
        for chrom in chrs:
            if chrom not in set(bed['chr']):
                print(f'Ignore {chrom}, since it does not appear in the BED file')
                continue
            profiles = process_chr(bamfile, bed, sample, chrom, sizes[chrom], fmin, fmax, window)
            results.append(profiles)

    return pd.concat(results, ignore_index=True)

if __name__ == '__main__':
    print('This is SCORE script')

    args = docopt(__doc__)

    sample = args['<sample_id>']
    fmin = int(args['--fmin'])
    fmax = int(args['--fmax'])
    input_dir = args['--input']
    output_dir = args['--output']

    bed = pd.read_table(bedpath, header=None)
    bed.columns = ['chr', 'start', 'end', 'gene', 'category']
    bed['chr'] = bed['chr'].str.replace(r'^Chr', 'chr', regex=True)
    bed['region_id'] = bed['chr'].astype(str) + ':' + \
                       bed['start'].astype(str) + ':' + \
                       bed['end'].astype(str) + ':' + \
                       bed['gene'].astype(str) + ':' + \
                       bed['category'].astype(str)
    chroms = bed['chr'].unique().tolist()

    results = main(sample, chroms, bed, fmin, fmax, window, input_dir, output_dir)

    output_path= os.path.join(output_dir, sample)
    os.makedirs(output_path, exist_ok=True)
    
    outfile = os.path.join(output_path, f"{sample}_profiles_updated.csv.gz")
    results.to_csv(outfile, compression='gzip', index=False)
