#!/usr/bin/env python3
"""
controltool.py

Count exact occurrences (including overlapping) of a set of motifs listed in a text file
across sequences in a FASTA file. This tool streams the FASTA, builds a single Aho–Corasick
automaton in memory from the motif list, and reports per-sequence or aggregate counts.

Usage examples:
  python controltool.py sequences.fasta motifs.txt -o counts.tsv
  python controltool.py seqs.fa motifs.txt --per-seq --revcomp --format json

Dependencies: Biopython (for SeqIO). No external Aho–Corasick dependency required.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import deque, defaultdict
from typing import Dict, Iterable, Iterator, List, Tuple

try:
    from Bio import SeqIO
except Exception as e:
    print("Biopython is required (pip install biopython)", file=sys.stderr)
    raise


# ---------------------- Aho-Corasick implementation (pure Python) ----------------------
class AhoNode:
    __slots__ = ("children", "fail", "outputs")

    def __init__(self):
        self.children: Dict[str, AhoNode] = {}
        self.fail: AhoNode | None = None
        self.outputs: List[str] = []


class AhoAutomaton:
    """Simple Aho–Corasick automaton. Patterns are strings.

    The automaton yields matches as tuples (pattern, start, end) where `end` is exclusive and
    `start` is inclusive. Overlapping matches are returned.
    """

    def __init__(self, patterns: Iterable[str]):
        self.root = AhoNode()
        self._build_trie(patterns)
        self._build_failures()

    def _build_trie(self, patterns: Iterable[str]):
        for pat in patterns:
            node = self.root
            for ch in pat:
                node = node.children.setdefault(ch, AhoNode())
            node.outputs.append(pat)

    def _build_failures(self):
        q = deque()
        # set fail for depth-1 nodes to root
        for ch, node in self.root.children.items():
            node.fail = self.root
            q.append(node)

        while q:
            current = q.popleft()
            for ch, child in current.children.items():
                q.append(child)
                # compute failure link
                f = current.fail
                while f is not None and ch not in f.children:
                    f = f.fail
                child.fail = f.children[ch] if (f and ch in f.children) else self.root
                child.outputs += child.fail.outputs if child.fail else []

    def iter_matches(self, text: str) -> Iterator[Tuple[str, int, int]]:
        node = self.root
        for i, ch in enumerate(text):
            while node is not None and ch not in node.children:
                node = node.fail
            if node is None:
                node = self.root
                continue
            node = node.children[ch]
            # outputs: list of patterns that end at position i
            for pat in node.outputs:
                start = i - len(pat) + 1
                yield pat, start, i + 1


# ---------------------- Utility functions ----------------------

def load_motifs(path: str) -> List[str]:
    """Load motifs from a text file: one motif per line, ignore blanks and comments (#).
    Returns motifs as uppercase strings, preserving input order but removing duplicates.
    """
    seen = set()
    motifs: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            motif = line.upper()
            if motif in seen:
                continue
            if len(motif) == 0:
                continue
            # Basic validation: only A/C/G/T allowed
            if any(ch not in "ACGT" for ch in motif):
                print(f"Warning: motif contains non-ACGT characters and will be skipped: {motif}", file=sys.stderr)
                continue
            motifs.append(motif)
            seen.add(motif)
    return motifs


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ---------------------- Processing and IO ----------------------

def process_fasta(
    fasta_path: str,
    automaton: AhoAutomaton,
    motif_to_label: Dict[str, str],
    per_seq: bool = False,
    min_len: int = 0,
    max_len: int = 10 ** 9,
    min_repeats: int = 1,
    count_units: bool = False,
) -> Tuple[Dict[str, int], Iterator[Tuple[str, str, int]]]:
    """Stream the FASTA file and use automaton to count matches.

    Returns:
      - totals: dict[label -> total_count]
      - if per_seq True, yields rows (label, seq_id, count) as they are completed for streaming output
    """
    totals: Dict[str, int] = defaultdict(int)

    def per_seq_generator() -> Iterator[Tuple[str, str, int]]:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            seq = str(rec.seq).upper()
            L = len(seq)
            if L < min_len or L > max_len:
                continue
            # collect start positions per label so we can detect consecutive repeats
            starts_by_label: Dict[str, List[int]] = defaultdict(list)
            for pat, start, end in automaton.iter_matches(seq):
                label = motif_to_label.get(pat, pat)
                starts_by_label[label].append(start)

            # analyze runs for each label
            for label, starts in starts_by_label.items():
                if not starts:
                    continue
                # starts are yielded in-order by the automaton; ensure sorted
                starts.sort()
                motif_len = None
                # find motif length: any pattern mapped to this label will have same length
                # find first pattern text with this label from motif_to_label (reverse lookup)
                # motif_to_label maps pattern->label, so find a pattern with this label
                for pat in motif_to_label:
                    if motif_to_label[pat] == label:
                        motif_len = len(pat)
                        break
                if motif_len is None:
                    # fallback: use diff between first two starts if possible
                    motif_len = starts[1] - starts[0] if len(starts) > 1 else 1

                run_count = 0
                run_len = 1
                prev = starts[0]
                for s in starts[1:]:
                    if s == prev + motif_len:
                        run_len += 1
                    else:
                        # finish previous run
                        if run_len >= min_repeats:
                            run_count += (run_len if count_units else 1)
                        run_len = 1
                    prev = s
                # last run
                if run_len >= min_repeats:
                    run_count += (run_len if count_units else 1)

                if run_count > 0:
                    totals[label] += run_count
                    yield label, rec.id, run_count

    if per_seq:
        return totals, per_seq_generator()

    # aggregate-only mode
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper()
        L = len(seq)
        if L < min_len or L > max_len:
            continue
        # collect starts per label for this record
        starts_by_label: Dict[str, List[int]] = defaultdict(list)
        for pat, start, end in automaton.iter_matches(seq):
            label = motif_to_label.get(pat, pat)
            starts_by_label[label].append(start)

        for label, starts in starts_by_label.items():
            if not starts:
                continue
            starts.sort()
            motif_len = None
            for pat in motif_to_label:
                if motif_to_label[pat] == label:
                    motif_len = len(pat)
                    break
            if motif_len is None:
                motif_len = starts[1] - starts[0] if len(starts) > 1 else 1

            run_count = 0
            run_len = 1
            prev = starts[0]
            for s in starts[1:]:
                if s == prev + motif_len:
                    run_len += 1
                else:
                    if run_len >= min_repeats:
                        run_count += (run_len if count_units else 1)
                    run_len = 1
                prev = s
            if run_len >= min_repeats:
                run_count += (run_len if count_units else 1)

            if run_count > 0:
                totals[label] += run_count

    return totals, iter(())


def write_output_tsv_aggregate(path: str, totals: Dict[str, int]):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("motif\ttotal_count\n")
        for motif in sorted(totals.keys()):
            fh.write(f"{motif}\t{totals[motif]}\n")


def write_output_tsv_perseq(path: str, rows: Iterable[Tuple[str, str, int]]):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("motif\tseq_id\tcount\n")
        for motif, seqid, cnt in rows:
            fh.write(f"{motif}\t{seqid}\t{cnt}\n")


def write_output_json(path: str, totals: Dict[str, int]):
    out = {"motifs": [{"motif": m, "total_count": totals[m]} for m in sorted(totals.keys())]}
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)


# ---------------------- CLI ----------------------

def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Count exact motif occurrences from a motif list in a FASTA.")
    parser.add_argument("fasta", help="path to input FASTA")
    parser.add_argument("motifs", help="path to motif text file (one motif per line)")
    parser.add_argument("-o", "--output", help="output path (default stdout)", default=None)
    parser.add_argument("--format", choices=["tsv", "json"], default="tsv", help="output format")
    parser.add_argument("--per-seq", action="store_true", help="emit counts per sequence (TSV only)")
    parser.add_argument("--revcomp", action="store_true", help="also count reverse-complements and add counts to the original motif")
    parser.add_argument("--min-length", type=int, default=0, help="skip sequences shorter than this")
    parser.add_argument("--max-length", type=int, default=10 ** 9, help="skip sequences longer than this")
    parser.add_argument("--min-repeats", type=int, default=1, help="minimum consecutive repeats to count as multiple occurrences")
    parser.add_argument("--count-units", action="store_true", help="count by units (e.g., each repeat) instead of runs")
    parser.add_argument("--stdout", action="store_true", help="write output to stdout instead of file (default if -o omitted)")
    args = parser.parse_args(argv)

    motifs = load_motifs(args.motifs)
    if not motifs:
        print("No valid motifs loaded. Exiting.", file=sys.stderr)
        return 2

    # Build pattern->label mapping. If revcomp, map both motif and its rc to the original motif label.
    motif_to_label: Dict[str, str] = {}
    for m in motifs:
        motif_to_label[m] = m
        if args.revcomp:
            rc = reverse_complement(m)
            motif_to_label[rc] = m

    # Build automaton with the set of keys
    automaton = AhoAutomaton(motif_to_label.keys())

    # Process FASTA
    totals, perseq_rows = process_fasta(
        args.fasta,
        automaton,
        motif_to_label,
        per_seq=args.per_seq,
        min_len=args.min_length,
        max_len=args.max_length,
        min_repeats=args.min_repeats,
        count_units=args.count_units,
    )

    outpath = args.output
    write_to_stdout = args.stdout or (outpath is None)

    if args.per_seq:
        # streaming write per-seq rows
        if args.format != "tsv":
            print("--per-seq currently supports only TSV output; using TSV.", file=sys.stderr)
        if write_to_stdout:
            # write header to stdout then rows
            sys.stdout.write("motif\tseq_id\tcount\n")
            for motif, seqid, cnt in perseq_rows:
                sys.stdout.write(f"{motif}\t{seqid}\t{cnt}\n")
        else:
            write_output_tsv_perseq(outpath, perseq_rows)
        # totals also available in 'totals' if needed
        return 0

    # aggregate mode
    if args.format == "tsv":
        if write_to_stdout:
            sys.stdout.write("motif\ttotal_count\n")
            for motif in sorted(totals.keys()):
                sys.stdout.write(f"{motif}\t{totals[motif]}\n")
        else:
            write_output_tsv_aggregate(outpath, totals)
        return 0

    # json
    if args.format == "json":
        if write_to_stdout:
            json.dump({"motifs": [{"motif": m, "total_count": totals[m]} for m in sorted(totals.keys())]}, sys.stdout, indent=2)
            sys.stdout.write("\n")
        else:
            write_output_json(outpath, totals)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
