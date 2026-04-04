#!/usr/bin/env python3
"""
Benchmark harness that uses `repeat_generator.py` to create FASTA inputs,
runs `DNASeqAna.py` on them while measuring execution time, and plots results.

This version uses the updated repeat_generator.generate_fasta API which accepts
keyword parameters like num_sequences, min_motif_size, max_motif_size,
num_distinct_motifs, min_repeats, and repeat_prob.

Usage:
  python3 benchmark_main.py --out-dir ./bench_out --runs 5

It requires matplotlib. Install with: pip install matplotlib
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path
import statistics
import json
import math

# local imports (generator lives in same folder)
from repeat_generator import generate_fasta

# Prefer non-interactive backend so script can run headless
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    plt = None


def run_single(fasta_path: Path, main_py: Path, extra_args: list[str]) -> float:
    start = time.perf_counter()
    # call DNASeqAna.py as a subprocess to isolate and measure full runtime
    cmd = [sys.executable, str(main_py), str(fasta_path)] + extra_args
    subprocess.check_call(cmd)
    end = time.perf_counter()
    return end - start


def frange_int(start: int, end: int, steps: int):
    if steps <= 1:
        yield start
        return
    for i in range(steps):
        # Use geometric spacing to cover wide ranges more evenly
        frac = i / (steps - 1)
        val = int(round(start * (end / start) ** frac))
        yield val


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("./bench_out"))
    parser.add_argument("--main-py", type=Path, default=Path("DNASeqAna.py"))
    # Backwards-compatible: original --seq-lengths interpreted as per-sequence length(s)
    parser.add_argument("--seq-lengths", type=int, nargs="*", default=[1000],
                        help="per-sequence length(s) for generated FASTA records (default: 1000)")
    # New: iterate over sequence-length range while keeping seq_count fixed
    parser.add_argument("--seq-length-start", type=int, default=None, help="start sequence length (e.g., 100)")
    parser.add_argument("--seq-length-end", type=int, default=None, help="end sequence length (e.g., 10000)")
    parser.add_argument("--seq-length-steps", type=int, default=None, help="how many steps between start and end (inclusive)")
    parser.add_argument("--num-seqs", type=int, default=1, help="number of sequences per FASTA file (deprecated if using range)")

    # New: iterate over number-of-sequences range
    parser.add_argument("--seq-count-start", type=int, default=None, help="start number of sequences (e.g., 10000)")
    parser.add_argument("--seq-count-end", type=int, default=None, help="end number of sequences (e.g., 5000000)")
    parser.add_argument("--seq-count-steps", type=int, default=None, help="how many steps between start and end (inclusive)")

    # repeat_generator parameters
    parser.add_argument("--min-motif-size", type=int, default=4, help="minimum motif size to sample")
    parser.add_argument("--max-motif-size", type=int, default=8, help="maximum motif size to sample")
    parser.add_argument("--distinct-motifs", type=int, default=5, help="number of distinct motifs to generate")
    parser.add_argument("--min-repeats", type=int, default=3, help="minimum repeat count when inserting tandem repeats")
    parser.add_argument("--repeat-prob", type=float, default=0.2, help="probability that a sequence contains a tandem repeat")

    parser.add_argument("--runs", type=int, default=1, help="how many repeats per input size")
    parser.add_argument("--extra-args", nargs=argparse.REMAINDER, default=[], help="extra args passed to DNASeqAna.py")
    parser.add_argument("--cleanup", action="store_true", help="delete generated FASTA files after processing to save disk space")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel processes to use (default: 4)")
    args = parser.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    # Pass workers argument to extra_args for DNASeqAna.py
    args.extra_args.extend(["--workers", str(args.workers)])

    # Determine sequence counts to test
    if args.seq_count_start is not None and args.seq_count_end is not None and args.seq_count_steps:
        seq_counts = list(dict.fromkeys(frange_int(args.seq_count_start, args.seq_count_end, args.seq_count_steps)))
    else:
        # Use --num-seqs and --seq-lengths behaviour: iterate over seq-lengths but fixed num_seqs
        seq_counts = [args.num_seqs]

    # Determine sequence lengths to test (range overrides explicit list)
    if args.seq_length_start is not None and args.seq_length_end is not None and args.seq_length_steps:
        seq_lengths = list(dict.fromkeys(frange_int(args.seq_length_start, args.seq_length_end, args.seq_length_steps)))
    else:
        seq_lengths = args.seq_lengths

    results = {}
    # Use per-sequence lengths list; we'll iterate counts x lengths
    for seq_length in seq_lengths:
        for seq_count in seq_counts:
            label = f"len{seq_length}_count{seq_count}"
            runtimes = []
            for r in range(args.runs):
                fasta_file = out / f"input_{label}_run_{r}.fasta"
                print(f"Generating {fasta_file} (len={seq_length}, count={seq_count}) ...")
                # Call the updated generate_fasta API
                generate_fasta(
                    filename=str(fasta_file),
                    num_sequences=seq_count,
                    seq_length=seq_length,
                    min_motif_size=args.min_motif_size,
                    max_motif_size=args.max_motif_size,
                    num_distinct_motifs=args.distinct_motifs,
                    min_repeats=args.min_repeats,
                    repeat_prob=args.repeat_prob,
                )

                print(f"Running DNASeqAna.py on {fasta_file} ...")
                try:
                    t = run_single(fasta_file, args.main_py, args.extra_args)
                except subprocess.CalledProcessError as e:
                    print("DNASeqAna.py failed with exit code", e.returncode)
                    t = float('nan')
                print(f"run {r} finished in {t:.2f}s")
                runtimes.append(t)

                # Optionally remove the generated FASTA to save disk
                if args.cleanup:
                    try:
                        fasta_file.unlink()
                        print(f"Deleted temporary FASTA {fasta_file}")
                    except Exception as e:
                        print(f"Warning: failed to delete {fasta_file}: {e}")

            # filter out NaNs if runs failed
            clean = [x for x in runtimes if not (x != x)]
            results[label] = {
                "seq_length": seq_length,
                "seq_count": seq_count,
                "runtimes": runtimes,
                "mean": statistics.mean(clean) if clean else float('nan'),
                "median": statistics.median(clean) if clean else float('nan')
            }

            # Save intermediate results so long runs can be resumed/inspected
            with open(out / "results.json", "w", encoding="utf-8") as f:
                json.dump(results, f, indent=2)

    # plot
    if plt is None:
        print("matplotlib not available; skipping plot")
        return

    # Decide whether sequence length or count was varied
    seq_lengths_set = sorted({v['seq_length'] for v in results.values()})
    seq_counts_set = sorted({v['seq_count'] for v in results.values()})

    if len(seq_lengths_set) > 1:
        # Sequence length varied: plot runtime vs sequence length (x-axis = sequence length).
        # Produce one plot per sequence-count (keeps comparisons simple when multiple counts are tested).
        for seq_count in seq_counts_set:
            xs = []
            ys = []
            for k, v in sorted(results.items(), key=lambda kv: kv[1]['seq_length']):
                if v['seq_count'] != seq_count:
                    continue
                xs.append(v['seq_length'])
                ys.append(v['mean'])
            if not xs:
                continue
            plt.figure()
            plt.plot(xs, ys, marker='o')
            plt.xlabel('sequence length')
            plt.ylabel('time (s)')
            plt.title(f'DNASeqAna.py runtime vs sequence length (seq_count={seq_count})')
            plt.grid(True)
            outpng = out / f'runtime_count{seq_count}_vs_length.png'
            plt.savefig(outpng)
            print('Saved plot to', outpng)
    else:
        # Sequence count varied (or neither varied): plot runtime vs number of sequences per seq_length
        for seq_length in seq_lengths_set:
            xs = []
            ys = []
            for k, v in sorted(results.items(), key=lambda kv: kv[1]['seq_count']):
                if v['seq_length'] != seq_length:
                    continue
                xs.append(v['seq_count'])
                ys.append(v['mean'])
            if not xs:
                continue
            plt.figure()
            plt.plot(xs, ys, marker='o')
            plt.xscale('log')
            plt.xlabel('number of sequences')
            plt.ylabel('time (s)')
            plt.title(f'DNASeqAna.py runtime vs number of sequences (seq_length={seq_length})')
            plt.grid(True)
            outpng = out / f'runtime_len{seq_length}.png'
            plt.savefig(outpng)
            print('Saved plot to', outpng)


if __name__ == '__main__':
    main()
