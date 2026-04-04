
import random
import textwrap
import optparse


def generate_random_dna(length):
    return "".join(random.choices("ACGT", k=length))

def generate_fasta(
    filename="synthetic_reads.fasta",
    num_sequences=100,
    seq_length=150,
    min_motif_size=4,
    max_motif_size=8,
    num_distinct_motifs=5,
    min_repeats=3,
    repeat_prob=0.2,
):
    # Generate a random pool of motifs to pick from:
    motif_pool = []
    for _ in range(num_distinct_motifs):
        motif_len = random.randint(min_motif_size, max_motif_size)
        motif_pool.append(generate_random_dna(motif_len))

    with open(filename, "w") as f:
        for i in range(num_sequences):
            sequence = list(generate_random_dna(seq_length))
            is_tr = False

            # Roll the dice if this sequence will contain a tandem repeat
            if random.random() < repeat_prob:
                is_tr = True
                motif = random.choice(motif_pool)
                repeat_count = min_repeats
                tr_sequence = motif * repeat_count
                tr_len = len(tr_sequence)

                if tr_len <= seq_length:
                    start_pos = random.randint(0, seq_length - tr_len)
                    sequence[start_pos : start_pos + tr_len] = list(tr_sequence)
                else:
                    is_tr = False

            # Write to file & Wrap sequence to 60 chars (standard FASTA convention)
            label = f"seq_{i+1}"
            full_seq = "".join(sequence)
            f.write(f">{label}\n")
            f.write("\n".join(textwrap.wrap(full_seq, 60)) + "\n")

    print("Motifs: ", motif_pool)
    print(f"Successfully generated {num_sequences} sequences to {filename}")


def main():
    parser = optparse.OptionParser()
    parser.add_option('-o', '--output', dest="output", action='store', type="str", default="", help="Name for output file (Required).")
    parser.add_option('-n', '--num-seqs', dest="num_seqs", action='store', default=500, type="int", help="Number of DNA sequences to generate (Default: 500).")
    parser.add_option('-l', '--length', dest="length", action='store', default=150, type="int", help="Uniform length of generated sequences (Default: 150).")
    parser.add_option('-s', '--min-size', dest="min_size", action='store', default=4, type="int", help="Minimum size of repeating motif to sample (Default: 4).")
    parser.add_option('-x', '--max-size', dest="max_size", action='store', default=8, type="int", help="Maximum size of repeating motif to sample (Default: 8).")
    parser.add_option('-m', '--motifs', dest="motifs", action='store', default=5, type="int", help="Number of random tandem repeat motifs to sample (Default: 5).")
    parser.add_option('-r', '--repeats', dest="repeats", action='store', default=5, type="int", help="Minimum number of repeats in a row (Default: 5).")
    parser.add_option('-p', '--probability', dest="prob", action='store', default=0.1, type="float", help="Base probability for a sequence to contain tandem repeats (Default: 0.1)")

    opts, args = parser.parse_args()

    if len(opts.output) == 0:
        print("We need at least an output file name... (-o, --output)")
    else:
        generate_fasta(
            filename=opts.output,
            num_sequences=opts.num_seqs,
            seq_length=opts.length,
            min_motif_size=opts.min_size,
            max_motif_size=opts.max_size,
            num_distinct_motifs=opts.motifs,
            min_repeats=opts.repeats,
            repeat_prob=opts.prob
        )


if __name__ == "__main__":
    main()
