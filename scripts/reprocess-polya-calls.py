#!/usr/bin/env python3
#
# Copyright (c) 2025 Seoul National University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

import pod5
import click
import pandas as pd
import numpy as np
import subprocess as sp
from scipy.ndimage import percentile_filter
from itertools import groupby

SIGNAL_LOAD_MAX_LENGTH = 30000
POLYA_MAXIMUM_DEVIATION = 25
POLYA_SCORE_POSITIVE = 1
POLYA_SCORE_NEGATIVE = -1
POLYA_LOWER_LIMIT_SIGMA = 3
POLYA_START_LIMIT_SIGMA = 2

def read_bam_pod5_pair(bamfile, pod5file, batch_size=1000):
    outbuf = []
    reads_in_batch = []

    with pod5.Reader(pod5file) as pod5reader, \
            sp.Popen(['samtools', 'view', '-F4', '-h', bamfile],
                     stdout=sp.PIPE, bufsize=1024*1024) as bamreader:

        for line in bamreader.stdout:
            if len(outbuf) >= batch_size:
                yield outbuf, load_pod5_signal(pod5reader, reads_in_batch)
                outbuf.clear()
                reads_in_batch.clear()

            if line.startswith(b'@'):
                outbuf.append(line)
                continue

            outbuf.append(line)
            read_id = line.split(b'\t', 1)[0].decode()
            reads_in_batch.append(read_id)

        if len(outbuf) > 0:
            yield outbuf, load_pod5_signal(pod5reader, reads_in_batch)

def load_pod5_signal(pod5reader, read_ids):
    signals = {}

    for read in pod5reader.reads(read_ids):
        signals[str(read.read_id)] = read.signal[:SIGNAL_LOAD_MAX_LENGTH]

    return signals

def detect_polya_stretches(signal, polya_start, polya_median, polya_std, window_size=60):
    lowerpct = percentile_filter(signal, size=window_size, percentile=25)
    lower_threshold = polya_median - POLYA_LOWER_LIMIT_SIGMA * polya_std

    is_polya = (lowerpct > lower_threshold)
    stretches = []
    for is_polya_i, group in groupby(enumerate(is_polya[polya_start:]), key=lambda x: x[1]):
        group = list(group)
        in_polya = bool(group[0][1])
        stretches.append([group[0][0], group[-1][0] + 1, in_polya])

    stretches = pd.DataFrame(stretches, columns=['start', 'end', 'is_polya'])
    stretches['length'] = stretches['end'] - stretches['start']
    stretches['score'] = stretches['length'] * (
        stretches['is_polya'] * POLYA_SCORE_POSITIVE + (1 - stretches['is_polya']) * POLYA_SCORE_NEGATIVE)
    stretches['score_balance'] = stretches['score'].cumsum()

    return stretches

def refine_polya_call(bamline, signal):
    bamfields = bamline.split(b'\t')
    extra_tags = [tag.split(b':', 1) for tag in bamfields[11:]]
    extra_tags = {tag[0].decode(): tag[1] for tag in extra_tags}

    if 'pt' not in extra_tags or 'pa' not in extra_tags:
        # No poly(A) call from dorado
        return bamline

    pt = int(extra_tags['pt'][2:].decode())
    pa = list(map(int, extra_tags['pa'][4:].decode().split(',')))

    if pa[0] < 0 or pa[1] < 0:
        # Poly(A) call is failed in dorado
        return bamline

    dorado_anchor = pa[0]
    dorado_pAstart = pa[1]
    dorado_pAend = pa[2]
    dorado_pAlen = pt
    dorado_speed = dorado_pAlen / (dorado_pAend - dorado_pAstart)

    # Find the median and std of the poly(A) tail region from dorado call
    prelim_polya_signal = signal[dorado_pAstart:dorado_pAend]
    prelim_polya_median = np.median(prelim_polya_signal)
    prelim_polya_filtered = prelim_polya_signal[np.where(
        np.abs(prelim_polya_signal - prelim_polya_median) < POLYA_MAXIMUM_DEVIATION)]
    prelim_polya_std = np.std(prelim_polya_filtered)

    polya_low = prelim_polya_median - prelim_polya_std * POLYA_START_LIMIT_SIGMA

    # Find the poly(A) start from the signal
    scan_region_left = 1000
    scan_polya_validation_window = 200
    for pos in np.where(signal[scan_region_left:dorado_pAend] >= prelim_polya_median)[0]:
        pos += scan_region_left
        if signal[pos:pos+scan_polya_validation_window].mean() > polya_low:
            refined_pAstart = pos
            break
    else:
        raise ValueError('No poly(A) tail found')

    refined_pAstart = min(dorado_pAstart, refined_pAstart)

    # Find the poly(A) end from the signal
    stretches = detect_polya_stretches(signal, refined_pAstart, prelim_polya_median, prelim_polya_std)
    terminal_block = stretches.loc[stretches['score_balance'].idxmax()]
    refined_pAend = terminal_block.end + refined_pAstart
    refined_pAlen = round(dorado_speed * (refined_pAend - refined_pAstart))

    index_pt = [i for i, tag in enumerate(bamfields[11:], 11) if tag.startswith(b'pt:')][0]
    index_pa = [i for i, tag in enumerate(bamfields[11:], 11) if tag.startswith(b'pa:')][0]

    bamfields[index_pt] = f'pt:i:{refined_pAlen}'.encode()
    bamfields[index_pa] = f'pa:B:i,{dorado_anchor},{refined_pAstart},{refined_pAend},-1,-1'.encode()

    return b'\t'.join(bamfields)

@click.command()
@click.argument('bamfile', type=click.Path(exists=True))
@click.argument('pod5file', type=click.Path(exists=True))
@click.argument('bamout', type=click.Path())
def main(bamfile, pod5file, bamout):
    with sp.Popen(['samtools', 'view', '-b', '-o', bamout], stdin=sp.PIPE) as bamwriter:
        bamoutf = bamwriter.stdin

        for bamlines, reads in read_bam_pod5_pair(bamfile, pod5file, 100):
            for bamline in bamlines:
                if bamline.startswith(b'@'):
                    bamoutf.write(bamline)
                    continue

                read_id = bamline.split(b'\t', 1)[0].decode()
                signal = reads[read_id]

                refined_bamline = refine_polya_call(bamline, signal)
                bamoutf.write(refined_bamline)
        
        bamoutf.flush()
        bamoutf.close()

if __name__ == '__main__':
    main()