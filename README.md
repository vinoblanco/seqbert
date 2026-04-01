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
- (optional) psutil — improves available-memory detection

Install dependencies:

    pip install -r requirements.txt

Or at minimum:

    pip install biopython
    pip install psutil   # optional, recommended

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
    reverse_comp text NOT NULL,
    UNIQUE (seq_number, motif)
  )

- After processing, the database is aggregated and exported to the tab-separated `output.txt` file (or the path you provide). The exported columns include motif, repeats, occurrences, proportion and reverse complement statistics.

## Key functions (implementation notes)

- `fasta_in_chunks(fasta_path, max_ram_mb=None, ram_fraction=0.15, avg_bytes_per_base=1.5) -> Iterator[List[SeqRecord]]`
  - Adaptive chunk generator that sizes chunks according to available memory. If `max_ram_mb` is set it overrides the adaptive calculation; otherwise the function attempts to detect available RAM (via `psutil` if installed, or `/proc/meminfo` on Linux) and uses `ram_fraction` of it as a budget.

- `equal_fasta_chunks(path, chunk_size=5000) -> Iterator[List[SeqRecord]]`
  - Simple fixed-size chunking using `itertools.islice`. This is used by default by `main.py`.

- `process_sequence(seq_str: str) -> dict`
  - Builds a dictionary keyed by all dinucleotide pairs (AA, AC, ..., TT) where values are `LinkedList` instances containing positions where that dinucleotide occurs.

- `calculate_motif_repeat_in_sequence(ll: LinkedList, seq_str: str, motive_size: int) -> Iterator[(start, period, occurrences)]`
  - Scans a linked list of positions for a dinucleotide key and yields detected runs where a motif repeats periodically. Returns start index, period (distance), and number of occurrences.

- `calculate_difference(ll, motive_size)`
  - Helper that groups adjacent position differences and yields counters and start positions used by higher-level routines.

- `statistical_repeats(base_tuple, seq_str: str, seq_id, min_repeats, max_repeats, motive_size)`
  - Aggregates findings for all dinucleotide keys for a single sequence, filters by `min_repeats`/`max_repeats`, computes canonical motifs and prepares rows for insertion into the database.

- `worker_process_chunk(records, min_repeats, max_repeats, motive_size)`
  - Runs in worker processes: takes lightweight records (tuples of id and sequence string), processes each sequence, and returns a list of rows matching the database schema for insertion.

- `write_repeats_to_txt(db, output_path='output.txt')`
  - Reads aggregated results from the temporary SQLite database and writes a summarized tab-separated file. The function performs aggregation (combining motifs with their reverse complements) before writing the final table.

## Memory / performance notes

- If `psutil` is installed, the script will use it to detect available memory and the optional adaptive `fasta_in_chunks` can then size chunks automatically. Without `psutil` the script falls back to reading `/proc/meminfo` on Linux or uses a conservative default.
- The multiprocessing approach uses `ProcessPoolExecutor(max_workers=4)` by default and limits the number of in-flight futures to avoid excessive memory use. You can change the number of workers by editing the `max_workers` argument in `main.py`.

## Limitations & TODOs

- The script currently focuses on dinucleotide keys as anchors for repeat detection — this is a heuristic and may miss repeats that are not captured by dinucleotide boundaries.
- Some internal functions assume the input sequence lengths and motif sizes are reasonable; extremely long sequences or very small motif sizes may need parameter tuning.
- The temporary database is stored in a NamedTemporaryFile — if you need a persistent DB, modify `main.py` to open a permanent sqlite file instead.

## Example

Given an input FASTA, run:

    python3 main.py input.fasta -m 4 -l 3 -u 10 -o output.txt

This will create `output.txt` (tab-separated) that lists motifs, their repeat counts and summary statistics.

## License

Unspecified. Add a `LICENSE` file if needed.
