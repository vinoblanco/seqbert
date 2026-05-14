# DNASeqAna

Small Python tool to detect periodic DNA repeats (motifs) from .fastq and .fasta files. The tool parses FASTA input in chunks, builds position maps for dinucleotide tuples, detects periodic distance patterns and writes discovered repeats into a temporary SQLite database and a .csv output file.

## Features

- Stream FASTA reading in fixed-size chunks
- Per-dinucleotide position indexing and list based analysis
- Detection of periodic repeats (motifs) and their properties: start, period, repeat count
- Canonical motif normalization and reverse complement handling
- Results written to a temporary SQLite database and exported to a csv file
- Multiprocessing worker pool (ProcessPoolExecutor) to process chunks in parallel

## Requirements

- Python 3.8+
- Biopython

Install dependencies:

    pip install -r requirements.txt

## Quick usage

Run the script on a FASTA file:

    python3 DNASeqAna.py path/to/input.fasta

CLI options:

- positional `fasta`: path to input FASTA file
- `-m`, `--min_motive_size` (int, default 4): minimum motif size to report
- `-M`, `--max_motive_size` (int, default 15): maximum motif size to report
- `-r`, `--min_repeats` (int, default 3): minimum number of repeat occurrences
- `-o`, `--output` (str, default `output.csv`): path to exported tab-separated output
- `-w`, `--workers` (int, default 4): number of worker processes for parallel processing
- `-c`, `--chunk_size` (int, default 5000): size per chunk


Example:

    python3 DNASeqAna.py input.fasta -m 4 -r 3 -o output.csv

## Behavior and output

- The script reads the input FASTA file in chunks. By default it uses a fixed-chunk helper (`equal_fasta_chunks`) which yields a fixed number of records per chunk.
- Each chunk is converted to a lightweight format (ID + sequence string) and submitted to a worker process.
- Workers run `worker_process_chunk`, which calls `find_denucleotides` (builds per-dinucleotide List position maps) and `detect_repeats` (extracts repeat candidates) and returns rows to insert into the central SQLite database.
- Database schema (created in a temporary file):

  CREATE TABLE repeats (
    seq_number text NOT NULL,
    motif text NOT NULL,
    period integer NOT NULL,
    repeat integer NOT NULL,
    reverse_comp text NOT NULL
  )

- After processing, the database is aggregated and exported to `output.csv` file (or the path you provide). The exported columns include motif, repeats, occurrences, proportion, the reverse complement statistics and the combined statistics.

## Memory / performance notes

- The multiprocessing approach uses `ProcessPoolExecutor` and 4 workers by default and limits the number of in-flight futures to avoid excessive memory use. You can change the number of workers by using the CLI-option `--workers`.

## Limitations & TODOs

- The script currently focuses on dinucleotide keys as anchors for repeat detection — this is a heuristic and may miss repeats that are not captured by dinucleotide boundaries.
- Some internal functions assume the input sequence lengths and motif sizes are reasonable; extremely long sequences or very small motif sizes may need parameter tuning.
- The temporary database is stored in a NamedTemporaryFile — if you need a persistent DB, modify `DNASeqAna.py` to open a permanent sqlite file instead.
- The script currently does not handle ambiguous bases (e.g. N) in the input sequences — these may need to be filtered or handled specially depending on your use case.
