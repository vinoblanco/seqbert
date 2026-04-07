# DNASeqAna

Small Python tool to detect periodic DNA repeats (motifs) from FASTA sequences. The tool parses FASTA input in chunks, builds position maps for dinucleotide tuples, detects periodic distance patterns and writes discovered repeats into a temporary SQLite database and a tab-separated output file.

## Features

- Stream FASTA reading in memory-sized chunks (adaptive) or fixed-size chunks
- Per-dinucleotide position indexing and linked-list based analysis
- Detection of periodic repeats (motifs) and their properties: start, period, repeat count
- Canonical motif normalization and reverse complement handling
- Results written to a temporary SQLite database and exported to a tab-separated text file
- Multiprocessing worker pool (ProcessPoolExecutor) to process chunks in parallel

## Requirements

- Python 3.8+
- Biopython
- psutil

Install dependencies:

    pip install -r requirements.txt

## Quick usage

Run the script on a FASTA file:

    python3 main.py path/to/input.fasta

CLI options:

- positional `fasta`: path to input FASTA file
- `-m`, `--motive_size` (int, default 4): minimum motif size to report
- `-l`, `--min_repeats` (int, default 3): minimum number of repeat occurrences
- `-u`, `--max_repeats` (int, default 10): maximum number of repeat occurrences
- `-o`, `--output` (str, default `output.txt`): path to exported tab-separated output
- `-w`, `--workers` (int, default 4): number of worker processes for parallel processing

Example:

    python3 main.py input.fasta -m 4 -l 3 -u 10 -o output.txt

The script prints progress messages, stores intermediate results in a temporary SQLite database, and writes a summarized result table to the specified output file when finished.

## Behavior and output

- The script reads the input FASTA file in chunks. By default it uses a fixed-chunk helper (`equal_fasta_chunks`) which yields a fixed number of records per chunk. There is also an adaptive reader (`fasta_in_chunks`) that splits the input according to available RAM.
- Each chunk is converted to a lightweight format (ID + sequence string) and submitted to a worker process.
- Workers run `worker_process_chunk`, which calls `process_sequence` (builds per-dinucleotide LinkedList position maps) and `statistical_repeats` (extracts repeat candidates) and returns rows to insert into the central SQLite database.
- Database schema (created in a temporary file):

  CREATE TABLE repeats (
    seq_number text NOT NULL,
    motif text NOT NULL,
    period integer NOT NULL,
    repeat integer NOT NULL,
    reverse_comp text NOT NULL
  )

- After processing, the database is aggregated and exported to `output.csv` file (or the path you provide). The exported columns include motif, repeats, occurrences, proportion and reverse complement statistics.

## Memory / performance notes

- The multiprocessing approach uses `ProcessPoolExecutor(max_workers=4)` by default and limits the number of in-flight futures to avoid excessive memory use. You can change the number of workers by editing the `max_workers` argument in `main.py`.

## Limitations & TODOs

- The script currently focuses on dinucleotide keys as anchors for repeat detection — this is a heuristic and may miss repeats that are not captured by dinucleotide boundaries.
- Some internal functions assume the input sequence lengths and motif sizes are reasonable; extremely long sequences or very small motif sizes may need parameter tuning.
- The temporary database is stored in a NamedTemporaryFile — if you need a persistent DB, modify `main.py` to open a permanent sqlite file instead.

## Example

Given an input FASTA, run:

    python3 main.py input.fasta -m 4 -l 3 -u 10 -o output.csv

This will create `output.csv` that lists motifs, their repeat counts and summary statistics.

## License

Unspecified. Add a `LICENSE` file if needed.
