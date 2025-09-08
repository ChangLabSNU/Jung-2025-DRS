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

import click
from Bio import Seq, SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
import numpy as np
import pysam
import itertools
import gzip

MINIMUM_SCORE = 5

def prepare_idrefs_for_global_alignment(allsequences, names, fivepmarker):
    trimmed_idrefs = {}
    for name in names:
        seq = str(allsequences[name].seq)
        idstart = seq.find(fivepmarker) + len(fivepmarker)
        trimmed_idrefs[name] = seq[idstart:]
    return trimmed_idrefs

def find_marker(aligner, fivepmarker, seq):
    markeraln = aligner.align(fivepmarker, seq)[0]
    return markeraln.coordinates[1][-1]

def find_polya_end(aln, seq):
    # Get signal moves
    moves = aln.get_tag('mv')
    headtrim = int(aln.get_tag('ts'))
    stride = moves[0]
    moves = np.array(moves[1:])

    # Get poly(A) end position in signal
    pa_tag = aln.get_tag('pa')
    if pa_tag is None or pa_tag[2] < 0:
        return None
    polya_end = pa_tag[2]
    if polya_end < headtrim:
        # poly(A) is removed completely.
        return None

    # Calculate the position of bases in signal
    posmap = len(seq) - np.cumsum(moves)
    polya_end_pos = (polya_end - headtrim) // stride
    return posmap[polya_end_pos]


def get_trimmed_seq(bamfile, fivepmarker):
    aligner = PairwiseAligner(mode='local', open_gap_score=-3, extend_gap_score=-1)

    identifiers = []

    for aln in pysam.AlignmentFile(bamfile, 'rb'):
        if aln.is_secondary or aln.is_supplementary:
            continue

        seq = aln.query_sequence

        markerend = find_marker(aligner, fivepmarker, seq)
        polya_end = find_polya_end(aln, seq)
        #print(seq[markerend-10:markerend].lower() + seq[markerend:polya_end] +
        #      seq[polya_end:polya_end+10].lower())
        if polya_end is not None:
            identifier_seq = seq[markerend:polya_end]
        else:
            identifier_seq = seq[markerend:]

        identifiers.append([aln.query_name, identifier_seq])

    return pd.DataFrame(identifiers, columns=['read_id', 'identifier_seq'])

def compare_to_idrefs(identifiers, idrefs, minimum_length=10):
    aligner = PairwiseAligner(mode='global', open_gap_score=-4, extend_gap_score=-1.5,
                              match_score=2, mismatch_score=-3)
    scores = {}

    for _, idseq in identifiers.iterrows():
        seq = idseq['identifier_seq']
        if len(seq) < minimum_length:
            continue

        scores[idseq['read_id']] = {
            refname: aligner.align(idref, seq)[0].score
            for refname, idref in idrefs.items()
        }

        #scores[idseq['read_id']]['length'] = len(seq)

    return pd.DataFrame.from_dict(scores, orient='index')

@click.command()
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--ref', type=click.Path(exists=True))
@click.option('--fivepmarker', type=str)
@click.option('--refname', type=str, multiple=True)
@click.option('--output', type=click.Path())
def main(bamfile, ref, fivepmarker, refname, output):
    allsequences = SeqIO.to_dict(SeqIO.parse(ref, 'fasta'))
    idrefs = prepare_idrefs_for_global_alignment(allsequences, refname, fivepmarker)
    identifiers = get_trimmed_seq(bamfile, fivepmarker)
    idmatch_scores = compare_to_idrefs(identifiers, idrefs)

    bestmatches = idmatch_scores.iloc[:, :len(refname)].idxmax(axis=1)
    bestmatch_score = idmatch_scores.iloc[:, :len(refname)].max(axis=1)
    bestmatches = pd.DataFrame({'bestmatch': bestmatches,
                                'bestmatch_score': bestmatch_score})
    bestmatches = bestmatches[bestmatches['bestmatch_score'] >= MINIMUM_SCORE]
    bestmatches.to_csv(output, header=False, sep='\t')

if __name__ == '__main__':
    main()