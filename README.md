# Supplementary Code for Jung *et al.* (2025)

[![DOI](https://zenodo.org/badge/1052605709.svg)](https://doi.org/10.5281/zenodo.17077444)

Nanopore direct RNA sequencing analysis pipeline for poly(A) tail profiling and RNA identification.

## Requirements

- Python 3.8+
- CUDA-capable GPU (for Dorado basecalling)
- Snakemake workflow manager

Install Python dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run the complete pipeline:
```bash
snakemake -j 8
```

Dry run to preview execution:
```bash
snakemake -n
```

## Input Data

Place raw sequencing data in:
- `rawdata/{run_name}/` - The raw output containing POD5 and FASTQ files

## Output

Processed results in:
- `demux/{run_name}/identifiers-{set_name}.txt` - Final demultiplexed RNA identifiers with poly(A) information

## Analysis Notebooks

Jupyter notebooks for data visualization:
- `notebooks/plot-polya-len-dists.ipynb` - Visualize poly(A) tail length distributions
- `notebooks/plot-signal-mixed-tailing.ipynb` - Analyze mixed tailing patterns in signal data

## Citation

Jung et al. (2025) details to be updated.

## Author

Hyeshik Chang <hyeshik@snu.ac.kr>

## License

MIT License - See source files for details