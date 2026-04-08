import itertools
import sqlite3
import tempfile
import os
import argparse
import csv
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from typing import Union, Dict, Optional, Any
from itertools import islice

from line_profiler import profile


class Node:
    def __init__(self, data):
        """
        Initializes a Node with data and a reference to the next node.
        :param data: The data to store in the node.
        """
        self.data = data
        self.next = None


class LinkedList:
    def __init__(self):
        """
        Initializes an empty linked list with head and tail set to None.
        """
        self.head = None
        self.tail = None

    def append(self, value):
        """
        Appends a new node with the given value to the end of the linked list.
        :param value: The value to append.
        """
        new_node = Node(value)

        if self.head is None:
            # List is empty, head and tail are the same node
            self.head = new_node
            self.tail = new_node
        else:
            # Link the current tail to the new node
            self.tail.next = new_node
            # Move the tail reference to the new node
            self.tail = new_node

    def iterate(self):
        """
        Iterates through the linked list and prints the data of each node.
        """
        node = self.head
        while node is not None:
            print(node.data)
            node = node.next

    def length(self):
        """
        Calculates the length of the linked list.
        :return: The number of nodes in the linked list.
        """
        count = 0
        node = self.head
        while node is not None:
            node = node.next
            count += 1
        return count

    def pairwise(self):
        """
        Generates pairs of consecutive nodes in the linked list.
        :yield: A tuple of two consecutive nodes.
        """
        node = self.head
        while node and node.next:
            yield node, node.next
            node = node.next


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
    In this method, a dictionary is created, which is populated with dinucleotide as keys and their occurrences as
    values.
    :param seq_str: the sequence for which the dict is created
    :return: the dictionary
    """
    dinucleotides: Dict[str, Optional[LinkedList]] = {
        "".join(kombi): None
        for kombi in itertools.product(["A", "C", "G", "T"], repeat=2)
    }

    for i in range(len(seq_str) - 1):
        key = str(seq_str[i : i + 2])
        if dinucleotides[key] is None:
            ll = LinkedList()
            dinucleotides[key] = ll
            ll.append(i)
        else:
            dinucleotides[key].append(i)

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
    for key, value in dinucleotides.items():
        if value is None or value.length() < min_repeats:
            continue
        for start, period, occurrences in search_motif(
            value, seq_str, motive_size, covered
        ):
            if min_repeats <= occurrences <= max_repeats:
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


def calculate_difference(ll: LinkedList, motive_size: int):
    """
    Evaluates the linked list, which contains all occurrences of a dinucleotide. Calculates whether the distance between
    two occurrences is periodic.
    :param ll: Linked list with the occurrences of the dinucleotide
    :param motive_size: minimum motif size
    :return: when the occurrence starts (start) and how often it occurs (count)
    """
    counter = {}
    start = None
    last = None

    for n, n1 in ll.pairwise():
        current = n1.data - n.data

        if last is None:
            last = current
            counter = {current: 1}
            start = n.data
            continue

        if current != last or current < motive_size:
            if counter:
                yield counter, start
            counter = {current: 1}
            start = n.data
            last = current
        else:
            counter[current] = counter.get(current, 0) + 1
            last = current

    if counter:
        yield counter, start


def search_motif(ll: LinkedList, seq_str: str, motive_size: int, covered: bytearray):
    """
    Finds motifs in the sequence based on the linked list with the occurrences of a dinucleotide
    Avoids areas that have already been marked as part of a found repeat.
    :param ll: Linked list with the occurrences of the dinucleotide
    :param seq_str: the sequence as a string
    :param motive_size: minimum motif size
    :param covered: bytearray, marks indices that have already been covered by a found repeat
    :return: when the occurrence starts (start), the period (period) and how often it occurs (occurrences)
    """
    run_period = None
    run_start = None
    run_occurrences = 1
    run_motif = None

    for n, n1 in ll.pairwise():
        # If one of the two positions has already been considered, skip it
        if covered[n.data] or covered[n1.data]:
            # If a run is already active, end it and start a new one
            yield from yield_run(covered, run_occurrences, run_period, run_start)
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        current = n1.data - n.data

        # Termination condition for motifs that are too short
        if current < motive_size:
            yield from yield_run(covered, run_occurrences, run_period, run_start)
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        # extract motifs for comparison
        motif_n = seq_str[n.data : n.data + current]
        motif_n1 = seq_str[n1.data : n1.data + current]

        # Termination condition for incomplete motifs at the end of the sequence
        if len(motif_n) != current or len(motif_n1) != current:
            yield from yield_run(covered, run_occurrences, run_period, run_start)
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        # If no run is active yet
        if run_period is None:
            # Start a new run if the motifs match
            if motif_n == motif_n1:
                run_period = current
                run_start = n.data
                run_occurrences = 2
                run_motif = motif_n
            # Reset values if no run can be started
            else:
                run_period = None
                run_start = None
                run_occurrences = 1
                run_motif = None

        # If a run is already active
        else:
            # Continue the run if the period and motif match
            if current == run_period and motif_n1 == run_motif:
                run_occurrences += 1

            # End the run and start a new run if the following motifs match
            else:
                # End the run
                yield from yield_run(covered, run_occurrences, run_period, run_start)
                # new run if motifs match
                if motif_n == motif_n1:
                    run_period = current
                    run_start = n.data
                    run_occurrences = 2
                    run_motif = motif_n
                else:
                    run_period = None
                    run_start = None
                    run_occurrences = 1
                    run_motif = None

    # Output the run from the last pass
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
    """
    Returns the reverse complement of a DNA sequence
    :param seq: DNA sequence
    :return: reverse complement of the DNA sequence
    """
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical_dna_motif(seq: str):
    """
    Sorts a sequence lexicographically and returns the smallest result
    :param seq: sequence
    :return: smallest sorting
    """
    candidates = [seq[i:] + seq[:i] for i in range(len(seq))]
    return min(candidates)


def worker_process_chunk(records, min_repeats, max_repeats, motive_size):
    """
    Prepares the data for processing in a process
    :param records: Sequence records only with id and sequence string
    :param min_repeats: minimum number of repetitions
    :param max_repeats: maximum number of repetitions
    :param motive_size: minimum motif size
    :return: data
    """
    rows = []
    for rec_id, seq_str in records:
        base_tuple = process_sequence(seq_str)
        repeats = statistical_repeats(
            base_tuple, seq_str, rec_id, min_repeats, max_repeats, motive_size
        )
        rows.extend(repeats)
    return rows


def write_output(
    db: Union[str, sqlite3.Connection], output_path: str = "output.txt"
) -> None:
    """
    Writes the found repeats from the database to a text file.
    :param output_path: Path to the output file.
    :param db: Database connection or path to the database file.
    """
    close_conn = False
    if isinstance(db, str):
        conn = sqlite3.connect(db)
        close_conn = True
    else:
        conn = db

    cur = conn.cursor()
    cur.execute("CREATE INDEX idx_motif ON repeats (motif)")
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

@profile
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

    starttime = time.time()
    with ProcessPoolExecutor(args.workers) as executor:
        print("Processing FASTA in chunks...")

        futures = set()

        for chunk in equal_fasta_chunks(args.fasta, args.chunk_size):
            lightweight = [(rec.id, str(rec.seq)) for rec in chunk]

            future = executor.submit(
                worker_process_chunk,
                lightweight,
                args.min_repeats,
                args.max_repeats,
                args.motive_size,
            )
            futures.add(future)

            if len(futures) >= MAX_IN_FLIGHT:
                done = next(as_completed(futures))
                futures.remove(done)

                rows = done.result()
                if rows:
                    cur.executemany(
                        "INSERT INTO repeats (seq_number, motif, period, repeat, reverse_comp) VALUES (?,?,?,?,?)",
                        rows,
                    )

        for future in as_completed(futures):
            rows = future.result()
            if rows:
                cur.executemany(
                    "INSERT INTO repeats (seq_number, motif, period, repeat, reverse_comp) VALUES (?,?,?,?,?)",
                    rows,
                )

        conn.commit()
    endtime = time.time()

    print(f"Writing results to {args.output} ... {endtime - starttime}")
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
