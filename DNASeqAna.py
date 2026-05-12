import sqlite3
import tempfile
import os
import argparse
import csv
import time
from collections import defaultdict
from compression import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
from Bio import SeqIO
from typing import Union, Any
from itertools import islice

DNA_TRANS_TABLE = str.maketrans("ACGT", "TGCA")


def equal_fasta_chunks(path, chunk_size):
    """
    Generator that yields chunks of FASTA records from the given file path.
    Each chunk contains a specified number of records, allowing for memory-efficient processing of large FASTA files.
    :param path: the path to the FASTA file
    :param chunk_size: the number of records in each chunk
    :return: a generator yielding lists of FASTA records
    """

    lower_path = path.lower()
    is_gz = lower_path.endswith(".gz")

    if lower_path.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")):
        fmt = "fasta"
    elif lower_path.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        fmt = "fastq"
    else:
        raise ValueError(f"Unsupported file type: {path}")

    opener = gzip.open if is_gz else open

    with opener(path, "rt") as handle:
        records = SeqIO.parse(handle, fmt)
        while True:
            chunk = list(islice(records, chunk_size))
            if not chunk:
                break
            yield chunk


def find_dinucleotides(seq_str: str):
    """
    Creates a dictionary populated with dinucleotides as keys and their occurrences as values.
    :param seq_str: the sequence as a string
    :return: a dictionary with dinucleotides as keys and their occurrences as values
    """
    dinucleotides = defaultdict(list)

    for i in range(len(seq_str) - 1):
        dinucleotides[seq_str[i : i + 2]].append(i)

    return dinucleotides


def detect_repeats(
    dinucleotides: dict,
    seq_str: str,
    seq_id: int,
    min_repeats: int,
    min_motive_size: int,
    max_motive_size: int,
):
    """
    Evaluates the dictionary and writes the found repeats to a database
    :param dinucleotides: the dictionary
    :param seq_str: the sequence of the dictionary as a string
    :param seq_id: the ID of the sequence
    :param min_repeats: minimum number of occurrences for it to be considered a repeat
    :param min_motive_size: minimum motif size
    :param max_motive_size: maximum motif size
    :return: an array with the found repeats and their properties (Seq_ID, motif, period, occurrences, reverse complement)
    """
    repeats = []
    covered = bytearray(len(seq_str))

    for key, positions in dinucleotides.items():
        if len(positions) < min_repeats:
            continue

        for start, period, occurrences in search_motif(
            positions, seq_str, min_motive_size, max_motive_size, covered
        ):
            if min_repeats <= occurrences:
                end = start + period * occurrences
                motif = str(canonical_dna_motif(seq_str[start : start + period]))
                if end <= len(seq_str):
                    repeats.append(
                        (
                            seq_id,
                            motif,
                            period,
                            occurrences,
                            canonical_dna_motif(reverse_complement(motif)),
                        )
                    )
    return repeats


def search_motif(positions: list, seq_str: str, min_motive_size: int, max_motive_size: int, covered: bytearray):
    """
    Searches for motifs in the given positions and yields the start, period and number of occurrences of each motif found.
    :param positions: the positions to search for motifs
    :param seq_str: the sequence as a string
    :param min_motive_size: minimum motif size
    :param max_motive_size: maximum motif size
    :param covered: a bytearray to mark positions that are already covered by found motifs
    :return: a generator yielding the start, period and number of occurrences of each motif found
    """
    run_period = None
    run_start = None
    run_occurrences = 1
    run_motif = None

    for n, n1 in zip(positions, positions[1:]):

        if covered[n] or covered[n1]:
            yield from yield_run(covered, run_occurrences, run_period, run_start)
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        current = n1 - n

        if current < min_motive_size or current > max_motive_size:
            yield from yield_run(covered, run_occurrences, run_period, run_start)
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        motif_n = seq_str[n : n + current]
        motif_n1 = seq_str[n1 : n1 + current]

        if len(motif_n) != current or len(motif_n1) != current:
            yield from yield_run(covered, run_occurrences, run_period, run_start)
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        if run_period is None:
            if motif_n == motif_n1:
                run_period = current
                run_start = n
                run_occurrences = 2
                run_motif = motif_n
            else:
                continue
        else:
            if current == run_period and motif_n1 == run_motif:
                run_occurrences += 1
            else:
                yield from yield_run(covered, run_occurrences, run_period, run_start)
                if motif_n == motif_n1:
                    run_period = current
                    run_start = n
                    run_occurrences = 2
                    run_motif = motif_n
                else:
                    run_period = None
                    run_start = None
                    run_occurrences = 1
                    run_motif = None

    yield from yield_run(covered, run_occurrences, run_period, run_start)


def yield_run(
    covered: bytearray,
    run_occurrences: int,
    run_period: Any | None,
    run_start: Any | None,
):
    """
    Yields the start, period and number of occurrences of a motif if it meets the criteria and marks the positions as covered.
    :param covered: a bytearray to mark positions that are already covered by found motifs
    :param run_occurrences: the number of occurrences of the current motif
    :param run_period: the period of the current motif
    :param run_start: the start position of the current motif
    :return: a generator yielding the start, period and number of occurrences of a motif if it meets the criteria
    """
    if run_period is not None and run_occurrences >= 2:
        start = run_start
        end = run_start + run_period * run_occurrences
        if not any(covered[i] for i in range(start, min(end, len(covered)))):
            for i in range(start, min(end, len(covered))):
                covered[i] = 1
            yield run_start, run_period, run_occurrences


def reverse_complement(seq: str):
    """
    Returns the reverse complement of a DNA sequence.
    :param seq: the DNA sequence as a string
    :return: the reverse complement of the DNA sequence as a string
    """
    return seq.translate(DNA_TRANS_TABLE)[::-1]


def canonical_dna_motif(seq: str):
    """
    Returns the canonical form of a DNA motif by finding the lexicographically smallest rotation of the sequence.
    :param seq: the DNA motif as a string
    :return: the canonical form of the DNA motif as a string
    """
    return min(seq[i:] + seq[:i] for i in range(len(seq)))


def worker_process_chunk(chunk, min_repeats, min_motive_size, max_motive_size):
    """
    Prepares the data for processing in a process
    :param chunk: Sequence records only with id and sequence string
    :param min_repeats: minimum number of repetitions
    :param min_motive_size: minimum motif size
    :param max_motive_size: maximum motif size
    :return: data
    """
    rows = []
    for rec_id, seq_str in chunk:
        base_dict = find_dinucleotides(seq_str)
        repeats = detect_repeats(
            base_dict, seq_str, rec_id, min_repeats, min_motive_size, max_motive_size
        )
        rows.extend(repeats)
    return rows


def write_output(
    db: Union[str, sqlite3.Connection], output_path: str = "output.csv"
) -> None:
    """
    Writes the results from the database to a CSV file, including calculations for proportions and combined repeats.
    :param db: the database connection or path to the database file
    :param output_path: the path to the output CSV file (default: "output.csv
    """
    close_conn = False
    if isinstance(db, str):
        conn = sqlite3.connect(db)
        close_conn = True
    else:
        conn = db

    cur = conn.cursor()

    print("Building indexes for final export...")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_motif ON repeats (motif)")
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT SUM(repeat) FROM repeats")
    grand_total = cur.fetchone()[0] or 1

    cur.execute("""
                CREATE
                TEMP TABLE agg_temp AS
                SELECT motif,
                       SUM(repeat)       AS total_repeats,
                       COUNT(seq_number) AS occurrences,
                       reverse_comp
                FROM repeats
                GROUP BY motif, reverse_comp
                """)

    cur.execute("CREATE INDEX idx_agg_motif ON agg_temp(motif)")
    conn.commit()

    cur = conn.cursor()
    query = """
            SELECT a.motif,
                   a.total_repeats,
                   a.occurrences,
                   ROUND(a.total_repeats * 1.0 / ?, 2)                                  AS proportion,
                   a.reverse_comp,
                   b.total_repeats                                                      AS rc_total_repeats,
                   b.occurrences                                                        AS rc_occurrences,
                   ROUND(b.total_repeats * 1.0 / ?, 2)                                  AS rc_proportion,
                   (a.total_repeats + COALESCE(b.total_repeats, 0))                     AS combined_repeats,
                   ROUND((a.total_repeats + COALESCE(b.total_repeats, 0)) * 1.0 / ?, 2) AS combined_proportion
            FROM agg_temp a
                     LEFT JOIN agg_temp b ON a.reverse_comp = b.motif
            WHERE a.reverse_comp IS NOT NULL
              AND (a.total_repeats > COALESCE(b.total_repeats, 0)
                OR (a.total_repeats = COALESCE(b.total_repeats, 0) AND a.motif < b.motif))
            ORDER BY combined_proportion DESC \
            """

    cur.execute(query, (grand_total, grand_total, grand_total))

    with open(output_path, "w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow([i[0] for i in cur.description])
        csv_writer.writerows(cur)

    if close_conn:
        conn.close()

def main():
    print("Starting processing...")

    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".db")
    os.close(tmp_fd)
    conn = sqlite3.connect(tmp_path)

    cur = conn.cursor()
    cur.execute("PRAGMA synchronous = OFF;")
    cur.execute("PRAGMA journal_mode = MEMORY;")
    cur.execute("PRAGMA temp_store = MEMORY;")
    cur.execute("""
                CREATE TABLE repeats
                (
                    seq_number   text    NOT NULL,
                    motif        text    NOT NULL,
                    period       integer NOT NULL,
                    repeat       integer NOT NULL,
                    reverse_comp text    NOT NULL
                )
                """)

    conn.commit()

    MAX_IN_FLIGHT = args.workers * 2
    chunk_counter = 0

    starttime = time.time()
    with ProcessPoolExecutor(args.workers) as executor:
        print(f"Processing FASTA in chunks with {args.workers} workers...")
        futures = set()

        for chunk in equal_fasta_chunks(args.fasta, args.chunk_size):

            lightweight_chunk = [(rec.id, str(rec.seq)) for rec in chunk]

            future = executor.submit(
                worker_process_chunk,
                lightweight_chunk,
                args.min_repeats,
                args.min_motive_size,
                args.max_motive_size,
            )
            futures.add(future)

            if len(futures) >= MAX_IN_FLIGHT:
                done, futures = wait(futures, return_when=FIRST_COMPLETED)
                for d in done:
                    rows = d.result()
                    if rows:
                        cur.executemany(
                            "INSERT INTO repeats (seq_number, motif, period, repeat, reverse_comp) VALUES (?,?,?,?,?)",
                            rows,
                        )
                        chunk_counter += 1

                if chunk_counter >= 20:
                    conn.commit()
                    chunk_counter = 0

        for future in as_completed(futures):
            rows = future.result()
            if rows:
                cur.executemany(
                    "INSERT INTO repeats (seq_number, motif, period, repeat, reverse_comp) VALUES (?,?,?,?,?)",
                    rows,
                )

        conn.commit()

    endtime = time.time()
    print(f"Extraction finished in {endtime - starttime:.2f}s")

    print(f"Writing results to {args.output}...")
    write_output(conn, args.output)

    conn.close()
    os.remove(tmp_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Finds and saves periodic repeats in DNA sequences from a FASTA file."
    )
    parser.add_argument("fasta", help="Path to the input FASTA file")
    parser.add_argument(
        "-m",
        "--min_motive_size",
        type=int,
        default=4,
        help="Minimum size of the motif (default: 4)",
    )
    parser.add_argument(
        "-M",
        "--max_motive_size",
        type=int,
        default=15,
        help="Maximum size of the motif (default: 15)",
    )
    parser.add_argument(
        "-r",
        "--min_repeats",
        type=int,
        default=3,
        help="Minimum number of repetitions (default: 3)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output.csv",
        help="Path to the output file (default: output.csv)",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=4,
        help="Number of parallel processes (default: 4)",
    )
    parser.add_argument(
        "-c",
        "--chunk_size",
        type=int,
        default=5000,
        help="Size of each chunk (default: 5000)",
    )
    args = parser.parse_args()

    main()
