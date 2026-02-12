import itertools
import sqlite3, tempfile
import os
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from linkedList import LinkedList
from typing import List, Iterator, Union
from itertools import islice

try:
    import psutil
except ImportError:
    psutil = None

#todo reverse complement berücksichtigen (evtl reverse comp schon mit in db schreiben)

#todo mit trfinder vergleichen
#todo stresstest mit 1 mo seq und hardware berücksichtigen

#todo einspeichern von motiven wenn sie zwei mal in einer sequenz vorkommen (auch wenn sie nicht periodisch sind)

def _get_available_memory_bytes() -> int:
    """
    Gibt den verfügbaren RAM in Bytes zurück. Bevorzuge psutil, wenn installiert,
    ansonsten lese /proc/meminfo (Linux).
    """
    if psutil:
        return int(psutil.virtual_memory().available)
    #Fallback for Linux
    meminfo_path = "/proc/meminfo"
    if os.path.exists(meminfo_path):
        with open(meminfo_path, "r") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    parts = line.split()
                    # value is in kB
                    return int(parts[1]) * 1024
    #Last resort: assume 256 MB available
    return 256 * 1024 * 1024

def fasta_in_chunks(fasta_path: str,
                    max_ram_mb: int = None,
                    ram_fraction: float = 0.5,
                    avg_bytes_per_base: float = 1.5) -> Iterator[List]:
    """
    Gibt Listen von SeqRecord Objekten zurück, die nach verfügbarem RAM dimensioniert sind.
    - fasta_path: Pfad zur FASTA Datei
    - max_ram_mb: maximaler RAM in MB (überschreibt ram_fraction, wenn gesetzt)
    - ram_fraction: Anteil des verfügbaren RAMs, der genutzt werden soll (0.0 - 1.0)
    - avg_bytes_per_base: geschätzte durchschnittliche Bytes pro Base in SeqRecord
    """
    if max_ram_mb is not None:
        allowed_bytes = int(max_ram_mb * 1024 * 1024)
    else:
        avail = _get_available_memory_bytes()
        allowed_bytes = int(avail * float(ram_fraction))

    # safety floor
    allowed_bytes = max(allowed_bytes, 1 * 1024 * 1024)  # at least 1 MB

    chunk = []
    chunk_bytes = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        rec_bases = len(record.seq)
        rec_estimated_bytes = int(rec_bases * avg_bytes_per_base)

        # If a single record is larger than allowed_bytes, yield it alone
        if rec_estimated_bytes >= allowed_bytes:
            if chunk:
                yield chunk
                chunk = []
                chunk_bytes = 0
            yield [record]
            continue

        # If adding this record would exceed the allowed budget, yield current chunk
        if chunk and (chunk_bytes + rec_estimated_bytes) > allowed_bytes:
            yield chunk
            chunk = []
            chunk_bytes = 0

        chunk.append(record)
        chunk_bytes += rec_estimated_bytes

    if chunk:
        yield chunk


def equal_fasta_chunks(path, chunk_size=5000):
    with open(path) as handle:
        records = SeqIO.parse(handle, "fasta")
        while True:
            chunk = list(islice(records, chunk_size))
            if not chunk:
                break
            yield chunk

def process_sequence(seq_str: str):
    """
    In dieser Methode wird ein Dictionary erstellt, welches mit Basentupeln als Keys und den Vorkommen dieser als
    Value befüllt wird.
    :param seq_str: die Sequenz, für welche das Dict erstellt wird
    :return: das Dictionary
    """
    base_tuple = {''.join(kombi): None for kombi in itertools.product(['A','C','G','T'], repeat=2)}
    #print(base_tuple)

    for i in range(len(seq_str) - 1):
        key = str(seq_str[i:i+2])
        if base_tuple[key] is None:
            ll = LinkedList()
            base_tuple[key] = ll
            ll.append(i)
        else:
            base_tuple[key].append(i)

    return base_tuple

#todo optimieren das nicht komplette Sequenz mehrfach analysiert wird
def statistical_repeats(base_tuple, seq_str: str, seq_id, min_repeats, max_repeats, motive_size):
    """
    Wertet das Dictionary aus und schreibt die gefundenen Repeats in eine Datenbank
    :param base_tuple: das Dictionary
    :param seq_str: die Sequenz des Dictionary als String
    :param seq_id: die ID der Sequenz
    :param min_repeats: mindest Anzahl der Vorkommen, damit es als Repeat angesehen wird
    :param max_repeats: maximale Anzahl der Vorkommen, damit es als Repeat angesehen wird
    :param motive_size: mindest Motivgröße
    :return: ein Array mit den gefundenen Repeats und deren Eigenschaften (Seq_ID, Motiv, Start, Periode, Anzahl)
    """
    repeats = []
    for key, value in base_tuple.items():
        #überspringt zu kurze oder leere Dictionary Einträge
        if value is None or value.lenght() < min_repeats:
            continue
        for start, period, occurrences in calculate_motif_repeat_in_sequence(value, seq_str, motive_size):
            if min_repeats <= occurrences <= max_repeats:
                end = start + period * occurrences
                motif = str(canonical_dna_motif(seq_str[start:start + period]))
                if end <= len(seq_str):
                    repeats.append((seq_id, motif, period, occurrences, canonical_dna_motif(reverse_complement(motif))))
    return repeats

def calculate_difference(ll, motive_size):
    """
    Wertet die Linked List aus, welche alle vorkommen eines Basentupels enthalten. Berechnet, ob der Abstand zwischen
    zwei Vorkommen periodisch ist.
    :param ll: Linked List mit den Vorkommen des Basentupels
    :param motive_size: mindest Motivgröße
    :return: wann das Vorkommen startet (start) und wie oft es auftritt (count)
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

def calculate_motif_repeat_in_sequence(ll, seq_str: str, motive_size: int):
    """
    Findet Motive in der Sequenz anhand der Linked List mit den Vorkommen eines Basentupels
    :param ll: Linked List mit den Vorkommen des Basentupels
    :param seq_str: die Sequenz als String
    :param motive_size: mindest Motivgröße
    :return: wann das Vorkommen startet (start), die Periode (period) und wie oft es auftritt (occurrences)
    """
    run_period = None
    run_start = None
    run_occurrences = 1
    run_motif = None

    for n, n1 in ll.pairwise():
        current = n1.data - n.data

        #Abbruchbedingung für zu kurze Motive
        if current < motive_size:
            if run_period is not None and run_occurrences >= 2:
                yield run_start, run_period, run_occurrences
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        #extrahiere Motive zum Vergleichen
        motif_n = seq_str[n.data:n.data + current]
        motif_n1 = seq_str[n1.data:n1.data + current]

        #Abbruchbedingung für unvollständige Motive am Ende der Sequenz
        if len(motif_n) != current or len(motif_n1) != current:
            if run_period is not None and run_occurrences >= 2:
                yield run_start, run_period, run_occurrences
            run_period = None
            run_start = None
            run_occurrences = 1
            run_motif = None
            continue

        #Wenn noch kein Run aktiv ist
        if run_period is None:
            #Starte neuen Run, wenn Motive übereinstimmen
            if motif_n == motif_n1:
                run_period = current
                run_start = n.data
                run_occurrences = 2
                run_motif = motif_n
            #Setze Werte zurück, wenn kein Run gestartet werden kann
            else:
                run_period = None
                run_start = None
                run_occurrences = 1
                run_motif = None

        #Wenn bereits ein Run aktiv ist
        else:
            #Run fortsetzen, wenn Periode und Motiv übereinstimmen
            if current == run_period and motif_n1 == run_motif:
                run_occurrences += 1

            #Run beenden und neuen Run starten, wenn die folgenden Motive übereinstimmt
            else:
                #Run beenden
                if run_occurrences >= 2:
                    yield run_start, run_period, run_occurrences
                #neuer Run, falls Motive übereinstimmen
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

    #Run vom letzten Durchlauf ausgeben
    if run_period is not None and run_occurrences >= 2:
        yield run_start, run_period, run_occurrences

def reverse_complement(seq):
    """
    Gibt das reverse complement einer DNA Sequenz zurück
    :param seq: DNA Sequenz
    :return: reverse complement der DNA Sequenz
    """
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def canonical_dna_motif(seq):
    """
    Sortiert eine Sequenz lexikographisch und gibt das kleinste Ergebnis zurück
    :param seq: Sequenz
    :return: kleinste Sortierung
    """
    candidates = [seq[i:] + seq[:i] for i in range(len(seq))]
    return min(candidates)

def worker_process_chunk(records, min_repeats, max_repeats, motive_size):
    """
    Vorbereitet die Daten für die Verarbeitung in einem Prozess
    :param records: Sequenz Datensätze nur mit id und Sequenz String
    :param min_repeats: minimale Anzahl der Wiederholungen
    :param max_repeats: maximale Anzahl der Wiederholungen
    :param motive_size: mindest Motivgröße
    :return: Daten
    """
    rows = []
    for rec_id, seq_str in records:
        base_tuple = process_sequence(seq_str)
        repeats = statistical_repeats(base_tuple, seq_str, rec_id , min_repeats, max_repeats, motive_size)
        rows.extend(repeats)
    return rows

def write_repeats_to_txt(db: Union[str, sqlite3.Connection], output_path: str = "output.txt") -> None:
    """
    Schreibt die gefundenen Repeats aus der Datenbank in eine Textdatei.
    :param output_path: Pfad zur Ausgabedatei.
    :param db: Datenbankverbindung oder Pfad zur Datenbankdatei.
    """
    close_conn = False
    if isinstance(db, str):
        conn = sqlite3.connect(db)
        close_conn = True
    else:
        conn = db

    cur = conn.cursor()
    #cur.execute("SELECT motif, " "SUM(repeat) AS total_repeats, " "COUNT(seq_number) AS occurrences, " "ROUND(SUM(repeat) * 1.0 / SUM(SUM(repeat)) OVER (), 2) AS proportion, " "reverse_comp " "FROM repeats GROUP BY motif ORDER BY total_repeats DESC")
    cur.execute("WITH aggregated AS ("
                "SELECT motif, "
                "SUM(repeat) AS total_repeats, "
                "COUNT(seq_number) AS occurrences, "
                "ROUND(SUM(repeat) * 1.0 / SUM(SUM(repeat)) OVER (), 2) AS proportion, "
                "reverse_comp "
                "FROM repeats "
                "GROUP BY motif, reverse_comp"
            "), paired AS ("
                "SELECT a.motif, "
                "a.total_repeats, "
                "a.occurrences, "
                "a.proportion, "
                "a.reverse_comp, "
                "b.motif AS rc_motif, "
                "b.total_repeats AS rc_total_repeats, "
                "b.occurrences AS rc_occurrences, "
                "b.proportion AS rc_proportion, "
                "(a.total_repeats + COALESCE(b.total_repeats, 0)) AS combined_repeats "
                "FROM aggregated a "
                "LEFT JOIN aggregated b ON a.reverse_comp = b.motif "
                "WHERE a.reverse_comp IS NOT NULL   AND (a.total_repeats > COALESCE(b.total_repeats, 0) OR (a.total_repeats = COALESCE(b.total_repeats, 0) AND a.motif < b.motif))"
            ") "
            "SELECT *, "
            "ROUND(combined_repeats * 1.0 / SUM(combined_repeats) OVER (), 2) "
            "AS combined_proportion "
            "FROM paired "
            "ORDER BY combined_proportion DESC")

    rows = cur.fetchall()

    with open(output_path, "w", encoding="utf-8") as f:
        # header
        f.write("motif\trepeats\toccurrences\tproportion\treverse_comp\trev_comp\trev_comp\trepeats\toccurrences\tproportion\ttotal_repeats\ttotal_proportions\n")
        for row in rows:
            f.write("\t".join(str(col) for col in row) + "\n")

    if close_conn:
        conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Findet und speichert periodische Repeats in DNA Sequenzen aus einer FASTA Datei.")
    parser.add_argument("fasta", help="Pfad zur Eingabe-FASTA Datei")
    parser.add_argument("-m", "--motive_size", type=int, default=4, help="Mindestgröße des Motivs (Standard: 4)")
    parser.add_argument("-l", "--min_repeats", type=int, default=3, help="Minimale Anzahl der Wiederholungen (Standard: 3)")
    parser.add_argument("-u", "--max_repeats", type=int, default=10, help="Maximale Anzahl der Wiederholungen (Standard: 10)")
    parser.add_argument("-o", "--output", type=str, default="output.txt", help="Pfad zur Ausgabedatei (Standard: output.txt)")
    args = parser.parse_args()

    print("Starting processing...")

    tmp = tempfile.NamedTemporaryFile(suffix=".db")
    conn = sqlite3.connect(tmp.name)

    cur = conn.cursor()
    cur.execute("""
    CREATE TABLE repeats (
        seq_number text NOT NULL,
        motif text NOT NULL,
        period integer NOT NULL,
        repeat integer NOT NULL,
        reverse_comp text NOT NULL,
        UNIQUE (seq_number, motif))
    """)

    cur.execute("CREATE INDEX idx_motif ON repeats (motif)")
    conn.commit()

    with ProcessPoolExecutor(max_workers=4) as executor:
        conn = sqlite3.connect(tmp.name)
        cur = conn.cursor()
        futures = []
        print("Processing FASTA in chunks...")
#        for chunk in equal_fasta_chunks(args.fasta):
        for chunk in fasta_in_chunks(args.fasta):
            lightweight = [(rec.id, str(rec.seq)) for rec in chunk]
            futures.append(executor.submit(worker_process_chunk, lightweight, args.min_repeats, args.max_repeats, args.motive_size))

            for future in as_completed(futures):
                rows = future.result()
                if rows:
                    cur.executemany("INSERT OR IGNORE INTO repeats VALUES (?,?,?,?,?)", rows)
                    conn.commit()

    print("Writing results to " + args.output + "...")
    write_repeats_to_txt(conn, args.output)
    conn.close()
    tmp.close()