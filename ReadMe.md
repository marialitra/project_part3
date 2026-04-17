
# Remote Homolog Search using Protein Embeddings & Approximate Nearest Neighbor Methods

## Authors

- **Lytra Maria** — [GitHub](https://github.com/marialitra)
- **Mylonaki Danai** — [GitHub](https://github.com/PoiLson)

## Main Idea

This project extends the work developed in previous assignments by applying Approximate Nearest Neighbor (ANN) methods to a real-world problem in computational biology.

In [Part 1](https://github.com/PoiLson/ann-search-algorithms-part1), we implemented classical ANN algorithms (LSH, Hypercube, IVF-Flat, IVFPQ), while in [Part 2](https://github.com/marialitra/ann-search-algorithms-part2) we introduced a learning-based approach using Neural LSH.

Building on these foundations, the current project tackles remote protein homology detection: identifying proteins with similar structure or function despite low sequence identity (<30% – the "Twilight Zone"), where tools like BLAST often fail.

We leverage protein embeddings from the pretrained ESM-2 model (facebook/esm2_t6_8M_UR50D) and adapt Approximate Nearest Neighbor (ANN) algorithms from prior assignments for efficient similarity search in embedding space.

Key features:

- **Embedding generation** from FASTA files via ESM-2 (mean pooling on last-layer representations).
- **ANN search** with methods: Euclidean LSH, Hypercube, IVF-Flat, IVFPQ, Neural LSH (or "all").
- **Evaluation** against BLAST baseline: Recall@N, QPS, and biological interpretation (e.g., remote homologs via UniProt/Pfam/GO annotations).
- **Hyperparameter tuning** via grid search for domain adaptation to protein embeddings.
- **Reports** with per-query summaries, neighbor tables, and bio comments.

The project is implemented in Python 3.10+ on Linux, integrating C binaries from Assignment 1 for ANN cores.

## Datasets

This repository already includes the protein datasets required for Part 3.

The dataset formats used by the reused ANN components from Parts 1 and 2 (MNIST and SIFT) are documented in [datasets.md](datasets.md). Please note that, while the MNIST files are included in the repository, the SIFT dataset must still be downloaded manually due to its size.

## Project Structure

- **`Algorithms/`**: ANN implementations from Assignments 1 & 2.
- **`Part3/`**: Protein-specific pipeline (ESM-2 embeddings, search benchmark, evaluation).

### Algorithms Directory

Contains shared ANN code reused for protein domain.

#### Subdirectories & Files in `Algorithms`

- **`AlgorithmsPart1/`**: C-based ANN cores (LSH, Hypercube, IVF-Flat, IVFPQ).
- **`Data/`**: Sample datasets (MNIST, SIFT) for testing prior algorithms.
- **`knngraphs/`**: KNN graphs output from IVFFlat.
- **`src/`**: Python wrappers and Neural LSH from Assignment 2.
- **`Makefile`**: Builds Python/C components for Assignment 2.

#### Source Files in `Algorithms/src/`

Modular Python scripts with documentation. Used for Neural LSH and wrappers.

- **`bruteforce.py`**: Runs brute-force KNN on train/query sets; saves results (.npy) and meta (params, avg time).
- **`check.py`**: Utility to verify versions of required libraries (for setup/dependency checks).
- **`datasetUtils.py`**: SIFT-specific functions for memmap-backed training.
- **`libraries.py`**: Imports common/custom libraries.
- **`neuralNet.py`**: Defines MLP/CNN classifiers; handles training for MNIST/SIFT.
- **`nlshBuild.py`**: Parses args, builds KNN graph, CSR, runs KaHIP, trains model, saves inverted index.
- **`nlshCore.py`**: Core query search and metrics computation for Neural LSH.
- **`nlshSearch.py`**: Parses args and invokes nlshCore functions.
- **`parseFiles.py`**: Reads IVFFlat output, MNIST/SIFT datasets.
- **`runSearchExe.py`**: Invokes C binary for IVFFlat via subprocess.
- **`utils.py`**: Parsing checks, filename generation, CSR building, training output saving/loading, data normalization.

#### Subdirectories & Files in `Algorithms/AlgorithmsPart1/`

C implementation for core ANN algorithms.

- **`Data/`**: MNIST/SIFT datasets.
- **`include/`**: Header files (see below).
- **`objectFiles/`**: Compiled .o files.
- **`src/`**: C source files (see below).
- **`Makefile`**: Builds the C executable.
- **`search`**: Compiled executable.
  - Used in Part 3 for ANN calls.
  - If it does not exist, it is automatically created by this project  
    (→ **no manual creation or changes are needed**)

#### Source Files in `Algorithms/AlgorithmsPart1/src/`

C implementations matching headers.

- **`datasets.c`**: Processes MNIST/SIFT datasets.
- **`hashtable.c`**: Fixed-capacity hash table with chaining.
- **`hypercube.c`**: Hypercube init, projection hashing, lookup, range search, destroy.
- **`ivfflat.c`**: IVFFlat clustering (KMeans++), assignment, lookup, range search, destroy.
- **`ivfpq.c`**: IVFPQ codebook training, residual computation, asymmetric distance lookup, range search, destroy.
- **`lsh.c`**: LSH init, amplified hashing, lookup, range search, destroy.
- **`main.c`**: Parses args, loads dataset, runs algorithm, frees memory.
- **`minheap.c`**: Min-heap ADT for top-N candidates.
- **`parseinput.c`**: Parses/validates command-line args.
- **`query.c`**: Executes queries, timings, brute-force cache, metrics, output saving.
- **`runAlgorithms.c`**: Initializes/runs each ANN algorithm.
- **`utils.c`**: Shared math (dot product, Euclidean dist, L2 norm, Hamming dist).

#### Header Files in `Algorithms/AlgorithmsPart1/include/`

Declarations for C functions/structures.

- **`datasets.h`**: Dataset utilities (MNIST/SIFT).
- **`hashtable.h`**: Hash table defs.
- **`hypercube.h`**: Hypercube structures/lookup.
- **`ivfflat.h`**: IVFFlat index/querying.
- **`ivfpq.h`**: IVFPQ encoding/search.
- **`lsh.h`**: LSH core defs.
- **`main.h`**: Library includes.
- **`minheap.h`**: Min-heap defs.
- **`parseinput.h`**: Arg parsing.
- **`query.h`**: Query execution/metrics.
- **`runAlgorithms.h`**: Algorithm runners.
- **`silhouette.h`**: Cluster evaluation (silhouette).
- **`utils.h`**: Common math ops.

### Part3 Directory

Protein homology pipeline.

#### Subdirectories & Files

- **`Data/`**: Protein datasets (swissprot.fasta, targets.fasta).
- **`output/`**: Generated files (see below).
- **`src/`**: Python scripts (see below).
- **`Makefile`**: Builds/runs Part3 components (e.g., BLAST).
- **`ReadMe.md`**: This file! (Project guide).
- **`requirements.txt`**: Required libraries (pip install -r).
- **`Biological Report.pdf`**: Report based on the biology of the project.
- **`Experimental Report.pdf`**: Report about the experiments throughout the project.

#### Output Subdirectories in `Part3/output/`

Organized results from pipeline steps.

- **`blast/`**: BLAST raw/processed outputs (e.g., tabular results, logs).
- **`embeddings/`**: Protein embeddings (.dat vectors, _ids.txt mappings).
- **`evaluation/`**: Markdown and helpers for biology report.
- **`search/`**: ANN results, final reports (e.g., per-query summaries).

#### Source Files in `Part3/src/`

ESM-2 integration and ANN adaptation for proteins.

- **`filter_blast.py`**: Filters BLAST (outfmt 6) by E-value/top-N (dedup subjects); ground truth TSV.
- **`generate_reports.py`**: Per-query/all-methods reports; Recall@N, QPS, tables with dist/identity/bio comments.
- **`libraries.py`**: Shared parsers, ESM loading, I/O, recall/QPS, report logic.
- **`parse_files.py`**: Arg parsers; ANN/BLAST result parsers to dicts.
- **`protein_embed.py`**: Generates ESM-2 embeddings from FASTA; .dat + _ids.txt.
- **`protein_search.py`**: Benchmark: query embeddings, ANN/BLAST runs, eval, reports.
- **`run_blast_methods.py`**: Runs/parses BLAST (via make); identity maps, QPS.
- **`run_methods.py`**: Orchestrates ANN (LSH etc.); embeddings, C calls, remapping.
- **`utils.py`**: BLAST filtering; ESM loading/embedding; ID remap; recall; bio comments.

## System Dependencies

This project requires NCBI BLAST+.

On Ubuntu / Debian-based systems, BLAST+ can be installed using:
```bash
sudo apt update
sudo apt install -y ncbi-blast+
```
Please ensure that the blastp, blastn, and related executables are available in your system PATH after installation.

## Installation Instructions

To install all required Python dependencies, run the following command inside the
Python environment you will use for this project:

```bash
pip install -r requirements.txt
```

## Usage

Use the `Part3/Makefile` for common workflows. It supports embedding generation on subsets, search benchmarks, BLAST runs, and grid search.

### Recommended Workflow (One Command)

```bash
# 1. Generate embeddings (only needed once, or when changing dataset)
make emb

# 2. Run full benchmark (ANN + BLAST + evaluation + reports)
make search
```
→ **`make search`** is the **main entry point** you should normally use.  
It calls `protein_search.py` with `-method all` by default and internally handles everything:

- Generates query embeddings (if not already present),
- Runs the selected ANN method(s),
- Automatically runs BLAST (if results are missing or outdated),
- Filters BLAST output to top-N **using exactly the same N** as specified in the search command,
- Computes Recall@N and QPS,
- Triggers detailed report generation.

This guarantees perfect consistency between ANN's top-N neighbors and BLAST's ground-truth top-N — you don't need to run `make blast` separately or worry about matching N values manually.

### Makefile Targets
- **`make emb`**: Generate embeddings for protein subset → output/embeddings/vectors.dat.
- **`make search`**: Run ANN search on vectors_50k.dat with queries (targets.fasta); method="all" → output/search/results.txt.
- **`make blast_db`**: Build BLAST database from swissprot_50k.fasta.
- **`make blast_search`**: Run blastp on targets.fasta → output/blast/search/blast_results.tsv.
- **`make blast`**: Full BLAST pipeline (db + search + filter top-N); outputs filtered TSV at output/blast/topN/blast_results_top${N}.tsv.
  - **Customize N**: Override via `make blast N=10` (default: 50). This changes top-N for filtering (used as ground truth for Recall@N).
  - E-value is hardcoded to 1e-2 for noise reduction; edit Makefile if needed.
- **`make grid_search`**: Run hyperparameter tuning; method="all" (override via `make grid_search method=lsh`).
  - Outputs per-config results in output/grid_search/.

### Direct Script Usage

- **Generate Embeddings**: `python Part3/src/protein_embed.py -i Data/swissprot.fasta -o output/embeddings/protein_vectors.dat`
- **Search Benchmark**: `python Part3/src/protein_search.py -d output/embeddings/protein_vectors.dat -q Data/targets.fasta -o output/search/results.txt -method <all|lsh|hypercube|neural|ivfflat|ivfpq>`
  - Note: We chose "neural" as the keyword for Neural LSH to align with the assignment's terminology.

### Output Behavior & Error Handling

When a method fails to complete:

- For single-method runs (`lsh`, `hypercube`, etc.): any existing `results.txt` is deleted.
- For `-method all`: failed methods are skipped — only successful ones contribute to the final report.
  - If nothing runs successfully: the output file is removed (no leftover empty file).

This prevents stale or misleading results from previous/failed executions.

### Customizing Parameters

- **Top-N (N) for ANN/Recall**: Use `-N <value>` flag in `protein_search.py` (default: 10 for neighbor tables; Recall@N uses Makefile's N for BLAST ground truth).
- **Algorithm Hyperparameters**: Defaults are tuned via grid search (best configs in code). Override via flags in `protein_search.py` (see parse_files.py for options):
  - LSH: `-k 2 -L 5 --lsh-w 20.0`
  - Hypercube: `-kproj 12 --hyper-w 20.0 --hyper-M 10000 -probes 100`
  - IVF-Flat: `--flat-kclusters 200 --flat-nprobe 100`
  - IVFPQ: `--pq-kclusters 500 --pq-nprobe 100 --pq-M 16 -nbits 8`
  - Neural LSH: `--nlsh-T 1000` (search); build params like `--nlsh-m 1800 --nlsh-layers 5` etc.
- Edit defaults directly in `Part3/src/parse_files.py` if needed (e.g., for persistent changes).

Outputs reports in `Part3/output/search/` with format as specified (e.g., [1] Method comparison table with Time/query, QPS, Recall@N; [2] Per-method neighbor tables with Rank, Neighbor ID, L2 Dist, BLAST Identity, In BLAST Top-N?, Bio comment).

## Version Control and Collaboration

The development of this project was managed using the Git version control system.

All source files, and experimental scripts were tracked through a dedicated Git repository to ensure collaborative development, change tracking, and reproducibility of results. The repository was hosted on a private GitHub project for version tracking and collaboration.
