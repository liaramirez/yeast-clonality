#!/usr/bin/env python3
"""
Figures and Table 1 for the Metschnikowia six-genome paper.

Run from /home/ldpramirez/nanopore_data (paths below are relative to it):

    python3 paper_figures.py --outdir figuras_paper \
        --pairwise-raw  results_raw/clonality/pairwise.tsv \
        --splits-raw    results_raw/clonality/group_splits.tsv \
        --pairwise-norm results_2kb/clonality/pairwise.tsv \
        --splits-norm   results_2kb/clonality/group_splits.tsv

Each figure is built independently: if an input is missing, that figure is
skipped with a message and the rest are still produced.
"""
import argparse, glob, os, re, subprocess, sys, statistics
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SAMPLES = ["Lev_1", "Lev_2", "Lev_3", "Lev_4", "Lev_5", "Lev_6"]
CLONE = ["Lev_1", "Lev_2", "Lev_3", "Lev_4", "Lev_5"]
CULTIVAR = {"Lev_1": "TA", "Lev_2": "GC", "Lev_3": "TA", "Lev_4": "ABA", "Lev_5": "PAIS", "Lev_6": "TA"}
BARCODE = {s: f"barcode0{i}" for i, s in enumerate(SAMPLES, 1)}
C_CLONE, C_OTHER, C_PUB, C_GREY = "#2b6c9e", "#c2572b", "#9a9a9a", "#d9d9d9"

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42})


def ok(path):
    return path and os.path.exists(path) and os.path.getsize(path) > 0


def skip(name, why):
    print(f"[omitida] {name}: {why}")


def save(fig, outdir, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"))
    plt.close(fig)
    print(f"[ok] {name}")


def clean(name):
    base = os.path.basename(name)
    m = re.search(r"GC[AF]_\d+\.\d+", base)
    if m:
        return m.group(0)
    return re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", base)


def read_fastani(path):
    d = {}
    for l in open(path):
        f = l.split()
        if len(f) < 3:
            continue
        try:
            d[(clean(f[0]), clean(f[1]))] = float(f[2])
        except ValueError:
            continue
    return d


# ---------------------------------------------------------------- Figure 1
def fig1(a, outdir):
    files = [a.ani_lev, a.ani_lev_refs, a.ani_refs]
    if not all(ok(f) for f in files):
        return skip("Fig1_ANI", "faltan tablas de fastANI (" + ", ".join(f for f in files if not ok(f)) + ")")
    d = {}
    for f in files:
        d.update(read_fastani(f))
    names = sorted({x for k in d for x in k})
    lev = [n for n in SAMPLES if n in names]
    pub = [n for n in names if n not in lev]
    order = lev + pub
    n = len(order)
    M = np.full((n, n), np.nan)
    for i, x in enumerate(order):
        for j, y in enumerate(order):
            if i == j:
                M[i, j] = 100
                continue
            v = [d[k] for k in ((x, y), (y, x)) if k in d]
            if v:
                M[i, j] = np.mean(v)
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(M, cmap="viridis", vmin=max(90, np.nanmin(M)), vmax=100)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(order, rotation=90, fontsize=6); ax.set_yticklabels(order, fontsize=6)
    ax.axhline(len(lev) - .5, color="w", lw=1); ax.axvline(len(lev) - .5, color="w", lw=1)
    fig.colorbar(im, ax=ax, label="ANI (%)", shrink=.8)
    ax.set_title("ANI among the six isolates and 20 public M. pulcherrima genomes")
    save(fig, outdir, "Fig1_ANI")

    pubpub = [M[i, j] for i in range(len(lev), n) for j in range(len(lev), n) if i < j and not np.isnan(M[i, j])]
    levpub = [M[i, j] for i in range(len(lev)) for j in range(len(lev), n) if not np.isnan(M[i, j])]
    with open(os.path.join(outdir, "Fig1_numbers.txt"), "w") as out:
        if pubpub:
            out.write(f"public vs public: median {np.median(pubpub):.2f}, range {min(pubpub):.2f}-{max(pubpub):.2f}\n")
        if levpub:
            best = [np.nanmax(M[i, len(lev):]) for i in range(len(lev))]
            out.write(f"isolates vs public, best hit per isolate: {min(best):.2f}-{max(best):.2f}\n")


# ---------------------------------------------------------------- Figure 2
def read_pairwise(path):
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for l in fh:
            r = dict(zip(header, l.rstrip("\n").split("\t")))
            key = tuple(sorted((r["sample_1"], r["sample_2"])))
            rows[key] = (float(r["discordance_pct"]), int(r["opposite_homozygotes"]))
    return rows


def read_splits(path):
    out = []
    with open(path) as fh:
        fh.readline()
        for l in fh:
            p, n = l.rstrip("\n").split("\t")
            out.append((p, int(n)))
    return out


def fig2(a, outdir):
    if not ok(a.ani_lev):
        return skip("Fig2", f"falta {a.ani_lev}")
    if not ok(a.pairwise_norm):
        return skip("Fig2", f"falta {a.pairwise_norm} (corre el pipeline con --chunk_size 2000)")
    ani = read_fastani(a.ani_lev)
    pw = read_pairwise(a.pairwise_norm)
    have_splits = ok(a.splits_raw) and ok(a.splits_norm)
    ncol = 3 if have_splits else 2
    fig, axes = plt.subplots(1, ncol, figsize=(3.4 * ncol, 3.2))

    # A: assembly ANI per pair
    pairs = sorted(pw)
    xs, ys, cs = [], [], []
    for p in pairs:
        v = [ani[k] for k in ((p[0], p[1]), (p[1], p[0])) if k in ani]
        if not v:
            continue
        within = p[0] in CLONE and p[1] in CLONE
        xs.append(np.mean(v)); ys.append(pw[p][1]); cs.append(C_CLONE if within else C_OTHER)
    ax = axes[0]
    for c, lab in ((C_CLONE, "within clone"), (C_OTHER, "clone vs Lev_6")):
        idx = [i for i, cc in enumerate(cs) if cc == c]
        ax.scatter([xs[i] for i in idx], [ys[i] + 1 for i in idx], c=c, s=28, label=lab, edgecolor="k", lw=.3)
    ax.set_yscale("log"); ax.set_xlabel("Assembly ANI (%)"); ax.set_ylabel("Opposite homozygotes + 1")
    ax.set_title("A  ANI vs genotypes", loc="left"); ax.legend(frameon=False, fontsize=7)

    # B: discordance per pair
    ax = axes[1]
    lab = [f"{p[0][-1]}-{p[1][-1]}" for p in pairs]
    disc = [pw[p][0] for p in pairs]
    col = [C_CLONE if (p[0] in CLONE and p[1] in CLONE) else C_OTHER for p in pairs]
    ax.bar(range(len(pairs)), disc, color=col)
    ax.set_yscale("log"); ax.set_xticks(range(len(pairs))); ax.set_xticklabels(lab, rotation=90, fontsize=6)
    ax.set_ylabel("Genotype discordance (%)"); ax.set_title("B  Pairwise discordance", loc="left")

    # C: splits raw vs normalized
    if have_splits:
        raw = dict(read_splits(a.splits_raw)); nor = dict(read_splits(a.splits_norm))
        keys = sorted(set(raw) | set(nor), key=lambda k: -raw.get(k, 0))
        ax = axes[2]; x = np.arange(len(keys)); w = .4
        ax.bar(x - w / 2, [raw.get(k, 0) for k in keys], w, color=C_GREY, label="original reads")
        ax.bar(x + w / 2, [nor.get(k, 0) for k in keys], w, color=C_CLONE, label="reads split to 2 kb")
        ax.set_xticks(x); ax.set_xticklabels([k.replace("Lev_", "") for k in keys], rotation=90, fontsize=6)
        ax.set_ylabel("Sites supporting the pair"); ax.set_title("C  Read-length artefact", loc="left")
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, outdir, "Fig2_clonality")


# ---------------------------------------------------------------- Figure 3
def af_values(vcf):
    cmd = ["bcftools", "query", "-i", 'GT="het" && FMT/DP>=10', "-f", "[%AD]\n", vcf]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    except Exception as e:
        return None
    v = []
    for l in out.splitlines():
        a = l.split(",")
        if len(a) >= 2 and a[0].isdigit() and a[1].isdigit():
            t = int(a[0]) + int(a[1])
            if t:
                v.append(int(a[1]) / t)
    return v


def fig3(a, outdir):
    vcfs = {s: f"ani/{s}.self.vcf.gz" for s in SAMPLES}
    if not any(ok(v) for v in vcfs.values()):
        return skip("Fig3", "faltan ani/Lev_*.self.vcf.gz")
    fig, axes = plt.subplots(2, 3, figsize=(8, 4.6), sharex=True)
    for ax, s in zip(axes.flat, SAMPLES):
        v = af_values(vcfs[s]) if ok(vcfs[s]) else None
        if not v:
            ax.set_title(f"{s}: no data"); continue
        ax.hist(v, bins=np.arange(0, 1.0001, .05), color=C_OTHER if s == "Lev_6" else C_CLONE)
        ax.set_title(f"{s} ({CULTIVAR[s]})", fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel("Alternative allele frequency")
    for ax in axes[:, 0]:
        ax.set_ylabel("Heterozygous sites")
    fig.tight_layout()
    save(fig, outdir, "Fig3_allele_frequencies")

    # optional: heterozygosity along the genome (needs ani/<s>.het.pos + ani/<s>.50kb.bed)
    if all(ok(f"ani/{s}.het.pos") and ok(f"ani/{s}.50kb.bed") for s in ("Lev_1", "Lev_6")):
        fig, axes = plt.subplots(2, 1, figsize=(7.5, 3.8), sharex=False)
        W = 50000
        for ax, s in zip(axes, ("Lev_1", "Lev_6")):
            L = {l.split()[0]: int(l.split()[2]) for l in open(f"ani/{s}.50kb.bed")}
            cnt = defaultdict(int)
            for l in open(f"ani/{s}.het.pos"):
                c, p = l.split()
                cnt[(c, int(p) // W)] += 1
            dens = []
            for c in sorted(L, key=lambda c: -L[c]):
                dens += [cnt[(c, w)] / W * 1000 for w in range(L[c] // W)]
            ax.bar(range(len(dens)), dens, width=1, color=C_OTHER if s == "Lev_6" else C_CLONE)
            ax.set_ylabel("het/kb"); ax.set_title(f"{s}: heterozygosity in 50-kb windows (contigs by size)", fontsize=8)
        axes[-1].set_xlabel("Window")
        fig.tight_layout()
        save(fig, outdir, "Fig3b_heterozygosity_profile")
    else:
        skip("Fig3b_heterozygosity_profile", "faltan ani/Lev_1.het.pos / Lev_6.het.pos (script del perfil de heterocigosis)")


# ---------------------------------------------------------------- Figure 4 + rDNA numbers
def asv2_pct(pileup, diag, good):
    n1 = n2 = 0
    for l in open(pileup):
        f = l.split("\t")
        pos = int(f[1])
        if pos not in good:
            continue
        ref = f[2].upper(); b1, b2 = diag[pos]
        s = re.sub(r"\^.", "", f[4]).replace("$", "")
        i = 0
        while i < len(s):
            c = s[i]
            if c in "+-":
                m = re.match(r"[+-](\d+)", s[i:]); i += len(m.group(0)) + int(m.group(1)); continue
            c2 = ref if c in ".," else c.upper()
            n1 += c2 == b1; n2 += c2 == b2; i += 1
    return (100 * n2 / (n1 + n2), (n1 + n2) / len(good)) if n1 + n2 else (None, 0)


def rdna(a):
    if not ok("ani/rdna/diag.tsv"):
        return {}
    diag = {int(l.split()[0]): tuple(l.split()[1:3]) for l in open("ani/rdna/diag.tsv")}
    good = set(a.rdna_positions)
    res = {}
    for s in SAMPLES:
        p = f"ani/rdna/{s}.pileup"
        if ok(p):
            res[s] = asv2_pct(p, diag, good)
    return res


def fig4(a, outdir, r):
    if not r:
        return skip("Fig4_rDNA", "faltan ani/rdna/diag.tsv o los pileups")
    fig, ax = plt.subplots(figsize=(3.6, 3))
    ss = [s for s in SAMPLES if s in r and r[s][0] is not None]
    ax.bar(range(len(ss)), [r[s][0] for s in ss], color=[C_OTHER if s == "Lev_6" else C_CLONE for s in ss])
    ax.set_xticks(range(len(ss))); ax.set_xticklabels([f"{s}\n{CULTIVAR[s]}" for s in ss], fontsize=7)
    ax.set_ylabel("rDNA copies carrying ASV2 (%)"); ax.set_ylim(0, 100)
    for i, s in enumerate(ss):
        ax.text(i, r[s][0] + 2, f"{r[s][0]:.1f}", ha="center", fontsize=7)
    save(fig, outdir, "Fig4_rDNA_ASV2")


# ---------------------------------------------------------------- Table 1 and mitochondria
def seqkit_n50(bc):
    fq = sorted(glob.glob(f"fastq_pass/{bc}/*.fastq.gz"))
    if not fq:
        return None, None
    try:
        out = subprocess.run(["seqkit", "stats", "-a", "-T", "-j", "4"] + fq, capture_output=True, text=True, check=True).stdout
    except Exception:
        return None, None
    lines = out.strip().splitlines()
    h = lines[0].split("\t"); tot = 0; n50s = []
    for l in lines[1:]:
        r = dict(zip(h, l.split("\t")))
        tot += int(r["sum_len"]); n50s.append((int(r["sum_len"]), float(r["N50"])))
    n50 = sum(a * b for a, b in n50s) / max(tot, 1)   # weighted approximation across files
    return tot, n50


def quast_table():
    for p in ("quast_out/report.tsv", "quast_out/transposed_report.tsv"):
        if ok(p):
            rows = [l.rstrip("\n").split("\t") for l in open(p)]
            if p.endswith("transposed_report.tsv"):
                h = rows[0]; return {clean(r[0]): dict(zip(h, r)) for r in rows[1:]}
            names = [clean(x) for x in rows[0][1:]]
            return {n: {r[0]: r[i + 1] for r in rows[1:]} for i, n in enumerate(names)}
    return {}


def busco_from_logs():
    res = {}
    for f in glob.glob("busco*.log") + glob.glob("busco_out/**/*.log", recursive=True) + glob.glob("busco_out/**/short_summary*.txt", recursive=True):
        txt = open(f, errors="ignore").read()
        for m in re.finditer(r"(Lev_\d)[^\n]*?\n?[^\n]*?C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\]", txt):
            res[m.group(1)] = (m.group(2), m.group(4))
    return res


def het_table():
    p = "ani/heterocigosis_resumen.tsv"
    if not ok(p):
        return {}
    rows = [l.rstrip("\n").split("\t") for l in open(p)]
    h = rows[0]
    return {r[0]: dict(zip(h, r)) for r in rows[1:]}


def table1(a, outdir, r):
    q = quast_table(); b = busco_from_logs(); het = het_table()
    cols = ["Isolate", "Cultivar", "Lineage", "Yield (Gb)", "Depth (x)", "Read N50 (kb)", "Assembly (Mb)", "GC (%)",
            "Contig N50 (kb)", "L50", "BUSCO C (%)", "BUSCO D (%)", "Het SNPs/kb", "ASV2 (%)"]
    lines = ["\t".join(cols)]
    for s in SAMPLES:
        tot, n50 = seqkit_n50(BARCODE[s]) if not a.skip_reads else (None, None)
        qq = q.get(s, {}) or q.get(f"{s}.medaka", {})
        def g(k):
            for kk in qq:
                if kk.lower().startswith(k.lower()):
                    return qq[kk]
            return ""
        size = g("Total length"); size = f"{int(size)/1e6:.2f}" if size.isdigit() else size
        n50c = g("N50"); n50c = f"{int(n50c)/1e3:.0f}" if n50c.isdigit() else n50c
        depth = f"{tot/float(g('Total length')):.0f}" if tot and g("Total length").isdigit() else ""
        bc = b.get(s, ("", ""))
        hk = het.get(s, {})
        hv = next((v for k, v in hk.items() if "kb" in k), "")
        asv = f"{r[s][0]:.1f}" if s in r and r[s][0] is not None else ""
        lines.append("\t".join([s, CULTIVAR[s], "clone" if s in CLONE else "second lineage",
                                f"{tot/1e9:.2f}" if tot else "", depth, f"{n50/1e3:.1f}" if n50 else "",
                                size, g("GC"), n50c, g("L50"), bc[0], bc[1], hv, asv]))
    open(os.path.join(outdir, "Table1.tsv"), "w").write("\n".join(lines) + "\n")
    print("[ok] Table1.tsv (revisa las celdas vacías: el script no encontró ese dato)")


def mito_table(outdir):
    fas = {s: f"ani/{s}_mt.fa" for s in SAMPLES}
    if not all(ok(f) for f in fas.values()):
        return skip("TableS2_mito", "faltan ani/Lev_*_mt.fa")
    lines = ["pair\tedit_distance\taligned_bp"]
    for i, x in enumerate(SAMPLES):
        for y in SAMPLES[i + 1:]:
            try:
                out = subprocess.run(["minimap2", "-cx", "asm5", fas[x], fas[y]], capture_output=True, text=True).stdout
            except FileNotFoundError:
                return skip("TableS2_mito", "minimap2 no está en el PATH (activa el ambiente)")
            nm = al = 0
            for l in out.splitlines():
                f = l.split("\t"); al += int(f[10])
                nm += sum(int(t.split(":")[2]) for t in f[12:] if t.startswith("NM:i:"))
            lines.append(f"{x}-{y}\t{nm}\t{al}")
    open(os.path.join(outdir, "TableS2_mito.tsv"), "w").write("\n".join(lines) + "\n")
    print("[ok] TableS2_mito.tsv")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="figuras_paper")
    p.add_argument("--ani-lev", default="ani/fastani.tsv")
    p.add_argument("--ani-lev-refs", default="ani/vs_refs.tsv")
    p.add_argument("--ani-refs", default="ani/refs_vs_refs.tsv")
    p.add_argument("--pairwise-norm", default="results_2kb/clonality/pairwise.tsv")
    p.add_argument("--pairwise-raw", default="results_raw/clonality/pairwise.tsv")
    p.add_argument("--splits-norm", default="results_2kb/clonality/group_splits.tsv")
    p.add_argument("--splits-raw", default="results_raw/clonality/group_splits.tsv")
    p.add_argument("--rdna-positions", type=int, nargs="+", default=[1578, 1585, 1586, 1587])
    p.add_argument("--skip-reads", action="store_true", help="no calcular rendimiento ni N50 de lecturas (más rápido)")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    fig1(a, a.outdir)
    fig2(a, a.outdir)
    fig3(a, a.outdir)
    r = rdna(a)
    fig4(a, a.outdir, r)
    with open(os.path.join(a.outdir, "rDNA_numbers.txt"), "w") as out:
        for s, (pct, depth) in r.items():
            out.write(f"{s}\tASV2 {pct:.1f}%\tmean reads per diagnostic site {depth:.0f}\n")
    table1(a, a.outdir, r)
    mito_table(a.outdir)
    print(f"\nListo. Revisa la carpeta {a.outdir}/")


if __name__ == "__main__":
    main()
