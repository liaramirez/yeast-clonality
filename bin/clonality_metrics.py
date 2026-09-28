#!/usr/bin/env python3
"""
Clonality metrics from a joint genotype table.

Input (--genotypes): tab-separated, first line = sample names, then one line per
SNP with one GT:DP field per sample (as written by `bcftools query -f '[%GT:%DP\t]'`).

For every pair of samples, on sites where both genotypes pass the depth filters:
  - discordance: % of sites where the two genotypes differ
  - opposite homozygotes: sites where one sample is 0/0 and the other 1/1
    (heterozygous-site mosaicism in assemblies or read-length effects cannot
    produce these; distinct strains produce thousands, clonal isolates ~0)

If --group lists >= 4 samples, it also counts, for every way of splitting the
group into 2 vs the rest, the sites where the 2 share one genotype and the rest
share a different one. Sub-lineages show as splits well above the others.

Outputs: pairwise.tsv, group_splits.tsv (if --group), summary.txt
"""
import argparse
import itertools
import statistics
import sys

GT_CODE = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genotypes", required=True)
    p.add_argument("--min-dp", type=int, default=10)
    p.add_argument("--max-dp-factor", type=float, default=2.0)
    p.add_argument("--group", default="", help="comma-separated sample names")
    return p.parse_args()


def read_table(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for i, r in enumerate(rows):
        if len(r) != len(header):
            sys.exit(f"Line {i + 2}: expected {len(header)} fields, got {len(r)}")
    return header, rows


def depth(field):
    d = field.split(":")[1] if ":" in field else "."
    return int(d) if d not in (".", "") else None


def main():
    a = parse_args()
    samples, rows = read_table(a.genotypes)
    n = len(samples)

    # median depth per sample, over sites where it has a depth value
    med = []
    for i in range(n):
        ds = [d for d in (depth(r[i]) for r in rows) if d is not None]
        med.append(statistics.median(ds) if ds else 0)

    def genotype(i, field):
        gt = field.split(":")[0]
        d = depth(field)
        if gt not in GT_CODE or d is None:
            return None
        if d < a.min_dp or d > a.max_dp_factor * med[i]:
            return None
        return GT_CODE[gt]

    G = [[genotype(i, f) for i, f in enumerate(r)] for r in rows]

    # pairwise
    pair = {}
    for x, y in itertools.combinations(range(n), 2):
        comp = disc = opp = 0
        for g in G:
            if g[x] is None or g[y] is None:
                continue
            comp += 1
            if g[x] != g[y]:
                disc += 1
                if {g[x], g[y]} == {0, 2}:
                    opp += 1
        pair[(x, y)] = (comp, disc, opp)

    with open("pairwise.tsv", "w") as out:
        out.write("sample_1\tsample_2\tcomparable_sites\tdiscordant\tdiscordance_pct\topposite_homozygotes\n")
        for (x, y), (c, d, o) in sorted(pair.items(), key=lambda kv: kv[1][1] / max(kv[1][0], 1)):
            out.write(f"{samples[x]}\t{samples[y]}\t{c}\t{d}\t{100 * d / max(c, 1):.3f}\t{o}\n")

    # group splits
    split_lines = []
    group = [s.strip() for s in a.group.split(",") if s.strip()]
    if group:
        missing = [s for s in group if s not in samples]
        if missing:
            sys.exit(f"--group samples not in the table: {', '.join(missing)}")
        idx = [samples.index(s) for s in group]
        if len(idx) >= 4:
            counts = {}
            for pair_idx in itertools.combinations(idx, 2):
                rest = [i for i in idx if i not in pair_idx]
                k = 0
                for g in G:
                    vals = [g[i] for i in idx]
                    if None in vals:
                        continue
                    A = [g[i] for i in pair_idx]
                    B = [g[i] for i in rest]
                    if A[0] == A[1] and len(set(B)) == 1 and A[0] != B[0]:
                        k += 1
                counts[pair_idx] = k
            with open("group_splits.tsv", "w") as out:
                out.write("pair\tsites\n")
                for p_, k in sorted(counts.items(), key=lambda kv: -kv[1]):
                    line = f"{samples[p_[0]]}+{samples[p_[1]]}\t{k}"
                    out.write(line + "\n")
                    split_lines.append(line)
            vals = sorted(counts.values())
            split_lines.append(f"(median {statistics.median(vals):.0f}, range {vals[0]}-{vals[-1]})")

    with open("summary.txt", "w") as out:
        out.write(f"Sites in table: {len(rows)}\n")
        out.write("Median depth: " + ", ".join(f"{s}={m:.0f}x" for s, m in zip(samples, med)) + "\n\n")
        out.write(f"{'pair':<28}{'comparable':>12}{'disc_%':>9}{'opp_hom':>9}\n")
        for (x, y), (c, d, o) in sorted(pair.items(), key=lambda kv: kv[1][1] / max(kv[1][0], 1)):
            out.write(f"{samples[x] + ' - ' + samples[y]:<28}{c:>12}{100 * d / max(c, 1):>9.2f}{o:>9}\n")
        if split_lines:
            out.write("\nGroup splits (2 vs rest):\n")
            out.write("\n".join("  " + l for l in split_lines) + "\n")
    print(open("summary.txt").read())


if __name__ == "__main__":
    main()
