import itertools
import sqlite3
import tempfile
import os
import argparse
import csv
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
from Bio import SeqIO
from typing import Union, Dict, Optional, Any
from itertools import islice

DNA_TRANS_TABLE = str.maketrans("ACGT", "TGCA")


def equal_fasta_chunks(path, chunk_size):
    with open(path) as handle:
        records = SeqIO.parse(handle, "fasta")
        while True:
            chunk = list(islice(records, chunk_size))
            if not chunk:
                break
            yield chunk


def process_sequence(seq_str: str):
    """
    Creates a dictionary populated with dinucleotides as keys and their occurrences as values.
    OPTIMIZATION: Uses built-in defaultdict and native C-arrays (lists) for instant appends.
    """
    dinucleotides = defaultdict(list)

    for i in range(len(seq_str) - 1):
        dinucleotides[seq_str[i: i + 2]].append(i)

    return dinucleotides


def statistical_repeats(
        dinucleotides: dict,
        seq_str: str,
        seq_id: int,
        min_repeats: int,
        max_repeats: int,
        motive_size: int,
):
    """
    Evaluates the dictionary and writes the found repeats to a database
    :param dinucleotides: the dictionary
    :param seq_str: the sequence of the dictionary as a string
    :param seq_id: the ID of the sequence
    :param min_repeats: minimum number of occurrences for it to be considered a repeat
    :param max_repeats: maximum number of occurrences for it to be considered a repeat
    :param motive_size: minimum motif size
    :return: an array with the found repeats and their properties (Seq_ID, motif, start, period, number)
    """
    repeats = []
    covered = bytearray(len(seq_str))

    for key, positions in dinucleotides.items():
        if len(positions) < min_repeats:
            continue

        for start, period, occurrences in search_motif(
                positions, seq_str, motive_size, covered
        ):
            if min_repeats <= occurrences <= max_repeats:
                end = start + period * occurrences
                motif = str(canonical_dna_motif(seq_str[start: start + period]))
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


def search_motif(positions: list, seq_str: str, motive_size: int, covered: bytearray):
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

        if current < motive_size:
            yield from yield_run(covered, run_occurrences, run_period, run_start)
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        motif_n = seq_str[n: n + current]
        motif_n1 = seq_str[n1: n1 + current]

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
                run_period = None
                run_start = None
                run_occurrences = 1
                run_motif = None
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
    if run_period is not None and run_occurrences >= 2:
        start = run_start
        end = run_start + run_period * run_occurrences
        if not any(covered[i] for i in range(start, min(end, len(covered)))):
            for i in range(start, min(end, len(covered))):
                covered[i] = 1
            yield run_start, run_period, run_occurrences


def reverse_complement(seq: str):
    return seq.translate(DNA_TRANS_TABLE)[::-1]


def canonical_dna_motif(seq: str):
    return min(seq[i:] + seq[:i] for i in range(len(seq)))


def worker_process_chunk(chunk, min_repeats, max_repeats, motive_size):
    """
    Prepares the data for processing in a process
    :param chunk: Sequence records only with id and sequence string
    :param min_repeats: minimum number of repetitions
    :param max_repeats: maximum number of repetitions
    :param motive_size: minimum motif size
    :return: data
    """
    rows = []
    for rec_id, seq_str in chunk:
        base_dict = process_sequence(seq_str)
        repeats = statistical_repeats(
            base_dict, seq_str, rec_id, min_repeats, max_repeats, motive_size
        )
        rows.extend(repeats)
    return rows


def write_output(
        db: Union[str, sqlite3.Connection], output_path: str = "output.csv"
) -> None:
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
                    period integer NOT NULL,
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
                args.max_repeats,
                args.motive_size,
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
        "--motive_size",
        type=int,
        default=4,
        help="Minimum size of the motif (default: 4)",
    )
    parser.add_argument(
        "-l",
        "--min_repeats",
        type=int,
        default=3,
        help="Minimum number of repetitions (default: 3)",
    )
    parser.add_argument(
        "-u",
        "--max_repeats",
        type=int,
        default=10,
        help="Maximum number of repetitions (default: 10)",
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
