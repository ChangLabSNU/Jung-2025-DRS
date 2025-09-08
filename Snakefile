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

import time
import os
import glob

RUNNAMES = ['rep1', 'rep2']

MINIMUM_READ_LENGTH = 200
REFERENCE_FASTA = 'references/ref.fa'

REPORTER_SEQS = ['A2', 'A7', 'control']
SPIKEIN_SEQS = ['QK01D', 'QK07D', 'QK08D', 'QK09D']
SPIKEIN_SUBSAMPLING = 10000
RNA_FIVEP_MARKER = 'GGCAAGTTGGACGCCCGCAAGATCCGCGAGATTCTCATTAAGGCCAAGAAGGGCGGCAAGATCGCCGTGTAA'

MODEL_HQ = 'rna004_130bps_sup@v5.2.0'
WDX_DEMUX_CONFIDENCE = 0.5

BARCODE_MAP = {
    'rep1': {
        'untreated': '4,7',
        'rg7834': '5,3',
    },
    'rep2': {
        'untreated': '4',
        'rg7834': '5',
    }
}
wildcard_constraints:
    run_name='[^_-]+',
    set_name='[^_-]+',
    sample='[^_-]+'

rule all:
    input:
        [f'demux/{run_name}/identifiers-{set_name}.txt'
         for run_name in RUNNAMES for set_name in ['reporter']]

rule collect_fastbasecall_sequences:
    output: 'fastbasecall/{run_name}/sequences.fastq.gz'
    threads: 16
    params:
        rundir='rawdata/{run_name}/'
    shell:
        '''
        find {params.rundir} -name '*.fastq.gz' -exec zcat {{}} \\; | \
        bgzip -c -@ {threads} > {output}
        '''

rule map_to_known_sequences:
    input:
        fastq='fastbasecall/{run_name}/sequences.fastq.gz',
        reference=REFERENCE_FASTA
    output:
        'alignments/{run_name}/primary-mapping.bam'
    threads: 32
    params:
        min_read_length=MINIMUM_READ_LENGTH,
        ref=REFERENCE_FASTA,
    shell:
        '''
        seqtk seq -L {params.min_read_length} {input.fastq} | \
        minimap2 -k 8 -w 4 -a -t {threads} --sam-hit-only {params.ref} - | \
        samtools view -bS -@ 4 -o {output}
        '''

rule extract_reporter_mapped_ids:
    input: 'alignments/{run_name}/primary-mapping.bam'
    output: 'alignments/{run_name}/reporter-mapped-ids.txt'
    run:
        seqpattern = '|'.join(REPORTER_SEQS)
        shell('samtools view -F 4 {input} | '
              'grep -E "\t({seqpattern})\t" | '
              'cut -f1 | sort -u > {output}')

rule extract_spikein_mapped_ids:
    input: 'alignments/{run_name}/primary-mapping.bam'
    output: 'alignments/{run_name}/spikein-mapped-ids.txt'
    run:
        seqpattern = '|'.join(SPIKEIN_SEQS)
        shell('samtools view -F 4 {input} | '
              'grep -E "\t({seqpattern})\t" | '
              'cut -f1 | sort -u | shuf | '
              'sed -n "1,{SPIKEIN_SUBSAMPLING}p" > {output}')

rule extract_pod5_subset:
    input: 'alignments/{run_name}/{set_name}-mapped-ids.txt'
    output: 'subset/{run_name}/{set_name}.pod5'
    threads: 8
    params:
        rundir='rawdata/{run_name}/pod5/'
    shell: 'pod5 filter --missing-ok -t {threads} -i {input} -o {output} {params.rundir}'

rule basecall_hq:
    input:
        pod5='subset/{run_name}/{set_name}.pod5',
        reference=REFERENCE_FASTA
    output: 'hqbasecall/{run_name}/{set_name}.bam'
    params:
        outputdir='hqbasecall/{run_name}/{set_name}'
    run:
        shell('rm -rf {params.outputdir}')
        shell('gpu -n 3 dorado basecaller -x cuda:all --reference {input.reference} \
              --mm2-opts "-k 8 -w 4" --estimate-poly-a -o {params.outputdir} \
              --emit-moves models/{MODEL_HQ} {input.pod5}')
        while True:
            try:
                shell('ls {params.outputdir}/*.bam')
            except Exception:
                print('Waiting for filesystem to sync...')
                time.sleep(5)
            else:
                break

        shell('ln {params.outputdir}/*.bam {output}')
        shell('samtools index {output}')

rule make_fastq:
    input: 'hqbasecall/{run_name}/{set_name}.bam'
    output: 'hqbasecall/{run_name}/{set_name}.fastq.gz'
    threads: 8
    shell: 'samtools fastq -F4 -F256 {input} | bgzip -c -@ {threads} > {output}'

rule refine_polya_calling:
    input:
        bam='hqbasecall/{run_name}/{set_name}.bam',
        pod5='subset/{run_name}/{set_name}.pod5'
    output: 'polya/{run_name}/{set_name}.bam'
    shell: 'python scripts/reprocess-polya-calls.py {input.bam} {input.pod5} {output}'

rule identify_rna_sequences:
    input: 'polya/{run_name}/{set_name}.bam'
    output: 'demux/{run_name}/rnaid-{set_name}.txt'
    shell: '''
        python scripts/detect-rna-identifiers.py --ref {REFERENCE_FASTA} \
        --fivepmarker {RNA_FIVEP_MARKER} --output {output} \
        --refname control --refname A2 --refname A7 {input}
    '''

rule run_barcode_demux:
    input: 'subset/{run_name}/{set_name}.pod5'
    output: 'demux/{run_name}/barcode-predictions-{set_name}.csv'
    params:
        outputdir='demux/{run_name}/{set_name}-warpdemux'
    threads: 8
    run:
        if os.path.exists(params.outputdir):
            shell('rm -rf ' + params.outputdir)

        shell('mamba run -n WDX warpdemux demux -i {input} -o {params.outputdir} \
                -m WDX4_rna004_v0_4_4 -j {threads}')

        outputfile = f'{params.outputdir}/warpdemux_*/predictions/barcode_predictions_0.csv'
        outputfiles = glob.glob(outputfile)
        assert len(outputfiles) == 1, 'Multiple output files found'

        shell('cat {outputfiles[0]} > {output}')

rule extract_read_ids_from_barcode_predictions:
    input: 'demux/{run_name}/barcode-predictions-{set_name}.csv'
    output: 'demux/{run_name}/read-ids-{set_name}-{sample}.txt'
    run:
        barcodes = BARCODE_MAP[wildcards.run_name][wildcards.sample].split(',')
        match_condition = ' || '.join([f'$2 == "{b}"' for b in barcodes])
        shell('awk -F, \'((' + match_condition + ') && $3 > {WDX_DEMUX_CONFIDENCE}) {{print $1 "\t" $3}}\' {input} > {output}')

def get_demux_read_id_files(wildcards):
    samples = list(BARCODE_MAP[wildcards.run_name].keys())
    return [f'demux/{wildcards.run_name}/read-ids-{wildcards.set_name}-{sample}.txt'
            for sample in samples]

rule merge_identifiers:
    input:
        rnaid='demux/{run_name}/rnaid-{set_name}.txt',
        barcodes=get_demux_read_id_files
    output: 'demux/{run_name}/identifiers-{set_name}.txt'
    run:
        import pandas as pd
        rnaid = pd.read_csv(input.rnaid, names=['read_id', 'rnaname', 'alnscore'],
                            sep='\t')
        barcodes = {}
        for barcodefile in input.barcodes:
            samplename = barcodefile.split('-')[-1].split('.')[0]
            df = pd.read_csv(barcodefile, names=['read_id', 'barcode_confidence'],
                             sep='\t')
            df['sample'] = samplename
            barcodes[samplename] = df
        barcodes = pd.concat(barcodes, axis=0).reset_index(drop=True)
        mgid = pd.merge(barcodes, rnaid, left_on='read_id', right_on='read_id')
        (mgid.sort_values(by=['sample', 'rnaname'])
             [['read_id', 'sample', 'rnaname', 'barcode_confidence', 'alnscore']]
             .to_csv(output[0], sep='\t', index=False))
