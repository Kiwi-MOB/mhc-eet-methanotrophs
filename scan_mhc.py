#!/usr/bin/env python3
"""
scan_mhc.py
-----------
Scan bacterial/archaeal proteomes for multiheme cytochrome c (MHC) candidates
by counting CXXCH / CXXCK heme-binding motifs.

Part of: "In silico survey of extracellular electron transfer potential in
aerobic methanotrophs and anammox bacteria"

USAGE
-----
    python scan_mhc.py --indir proteomes/ --outdir results/ --min-hemes 3

INPUT
-----
  proteomes/          a folder of protein FASTA files (.faa / .fasta / .fa),
                      ONE FILE PER GENOME.
                      Name each file with its assembly accession, e.g.
                        GCF_000008325.1.faa
                      (the filename stem becomes the genome_id)

  genomes.csv         OPTIONAL metadata table in the same folder as this
                      script. If present, its columns are joined onto the
                      output so your tables carry taxonomy. Must contain a
                      'genome_id' column matching the filename stems.

OUTPUT (written to --outdir)
----------------------------
  per_genome_summary.csv   one row per genome: proteome size, how many
                           candidates, max heme count. This is your Figure 1.
  candidates.csv           one row per candidate protein. This is your master
                           table -- you will add localization / domain /
                           neighbourhood columns to it by hand later.
  candidates.faa           FASTA of the candidate proteins ONLY.
                           <-- paste this into SignalP / DeepTMHMM /
                               InterProScan. Do NOT submit whole proteomes.

NOTES ON THE MOTIF SEARCH
-------------------------
  * CXXCH is the canonical c-type heme attachment motif; CXXCK is a known
    variant (e.g. in some hydroxylamine oxidoreductase-family proteins).
  * Matches are found with a LOOKAHEAD regex so that overlapping motifs
    (CXXCHXXCH) are both counted. A naive re.findall would miss these.
  * Chance occurrence is real: in a ~3000-protein proteome you will get some
    single-motif hits that are not cytochromes at all. That is exactly why
    the default threshold is >=3 and why localization filtering comes next.
  * This script does NOT prove anything is a cytochrome. It produces a
    candidate list. Treat every output row as a hypothesis.
"""

import argparse
import csv
import os
import re
import sys

# CXXCH / CXXCK, allowing overlaps via lookahead.
MOTIF_CXXCH = re.compile(r"(?=(C..CH))")
MOTIF_CXXCK = re.compile(r"(?=(C..CK))")

FASTA_EXT = (".faa", ".fasta", ".fa", ".pep")


def read_fasta(path):
    """Minimal FASTA parser. Yields (header, sequence). No dependencies."""
    header, chunks = None, []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def parse_header(header):
    """Split a FASTA header into (accession, description)."""
    parts = header.split(None, 1)
    acc = parts[0] if parts else header
    desc = parts[1] if len(parts) > 1 else ""
    return acc, desc


def load_metadata(path):
    """Load optional genomes.csv keyed on genome_id."""
    if not path or not os.path.exists(path):
        return {}, []
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return {}, []
    if "genome_id" not in rows[0]:
        sys.stderr.write(
            "WARNING: genomes.csv has no 'genome_id' column; ignoring it.\n")
        return {}, []
    extra_cols = [c for c in rows[0].keys() if c != "genome_id"]
    return {r["genome_id"]: r for r in rows}, extra_cols


def scan_proteome(path, min_hemes):
    """Scan one proteome file. Returns (summary_dict, list_of_candidates)."""
    genome_id = os.path.splitext(os.path.basename(path))[0]
    n_proteins = 0
    candidates = []
    max_hemes = 0

    for header, seq in read_fasta(path):
        n_proteins += 1
        seq = seq.upper().replace("*", "")

        n_cxxch = len(MOTIF_CXXCH.findall(seq))
        n_cxxck = len(MOTIF_CXXCK.findall(seq))
        n_total = n_cxxch + n_cxxck
        if n_total > max_hemes:
            max_hemes = n_total

        if n_total >= min_hemes:
            acc, desc = parse_header(header)
            candidates.append({
                "genome_id": genome_id,
                "protein_acc": acc,
                "description": desc,
                "length_aa": len(seq),
                "n_CXXCH": n_cxxch,
                "n_CXXCK": n_cxxck,
                "n_heme_motifs": n_total,
                "n_Cys": seq.count("C"),
                "sequence": seq,
            })

    summary = {
        "genome_id": genome_id,
        "n_proteins": n_proteins,
        "n_candidates": len(candidates),
        "max_heme_motifs": max_hemes,
        "candidates_per_1000_proteins": (
            round(1000 * len(candidates) / n_proteins, 2) if n_proteins else 0),
    }
    return summary, candidates


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True,
                    help="folder of protein FASTA files, one per genome")
    ap.add_argument("--outdir", default="results",
                    help="output folder (created if missing)")
    ap.add_argument("--min-hemes", type=int, default=3,
                    help="minimum CXXCH+CXXCK count to call a candidate (default 3)")
    ap.add_argument("--metadata", default="genomes.csv",
                    help="optional metadata CSV keyed on genome_id")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    files = sorted(
        os.path.join(args.indir, f)
        for f in os.listdir(args.indir)
        if f.lower().endswith(FASTA_EXT)
    )
    if not files:
        sys.exit(f"ERROR: no FASTA files ({', '.join(FASTA_EXT)}) found in {args.indir}")

    meta, extra_cols = load_metadata(args.metadata)

    summaries, all_candidates = [], []
    for path in files:
        summary, cands = scan_proteome(path, args.min_hemes)
        summaries.append(summary)
        all_candidates.extend(cands)
        print(f"{summary['genome_id']:<28} "
              f"{summary['n_proteins']:>6} proteins  "
              f"{summary['n_candidates']:>4} candidates  "
              f"max {summary['max_heme_motifs']:>3} motifs")

    # --- per-genome summary ---
    sum_cols = ["genome_id"] + extra_cols + [
        "n_proteins", "n_candidates", "candidates_per_1000_proteins",
        "max_heme_motifs"]
    sum_path = os.path.join(args.outdir, "per_genome_summary.csv")
    with open(sum_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sum_cols, extrasaction="ignore")
        w.writeheader()
        for s in summaries:
            row = dict(s)
            row.update({c: meta.get(s["genome_id"], {}).get(c, "") for c in extra_cols})
            w.writerow(row)

    # --- candidate table (sequence dropped; it lives in the .faa) ---
    # Blank columns are placeholders you fill in from the web tools later.
    manual_cols = ["signal_peptide", "predicted_localization",
                   "interpro_domains", "neighbour_porin", "notes"]
    cand_cols = ["genome_id"] + extra_cols + [
        "protein_acc", "description", "length_aa",
        "n_CXXCH", "n_CXXCK", "n_heme_motifs", "n_Cys"] + manual_cols
    cand_path = os.path.join(args.outdir, "candidates.csv")
    all_candidates.sort(key=lambda c: (-c["n_heme_motifs"], c["genome_id"]))
    with open(cand_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cand_cols, extrasaction="ignore")
        w.writeheader()
        for c in all_candidates:
            row = dict(c)
            row.update({col: meta.get(c["genome_id"], {}).get(col, "") for col in extra_cols})
            row.update({col: "" for col in manual_cols})
            w.writerow(row)

    # --- candidate FASTA, ready for SignalP / DeepTMHMM / InterProScan ---
    faa_path = os.path.join(args.outdir, "candidates.faa")
    with open(faa_path, "w") as fh:
        for c in all_candidates:
            fh.write(f">{c['genome_id']}|{c['protein_acc']} "
                     f"hemes={c['n_heme_motifs']} len={c['length_aa']}\n")
            seq = c["sequence"]
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    print("\n" + "=" * 60)
    print(f"genomes scanned : {len(summaries)}")
    print(f"candidates      : {len(all_candidates)} (>={args.min_hemes} heme motifs)")
    print(f"\nwrote {sum_path}")
    print(f"wrote {cand_path}")
    print(f"wrote {faa_path}")
    print("\nNext: submit candidates.faa to SignalP 6.0 and DeepTMHMM,")
    print("then fill the blank columns in candidates.csv.")


if __name__ == "__main__":
    main()
