#!/usr/bin/env python3
"""
merge_psortb.py
---------------
Parse PSORTb 3.0 "Long" format output and join it onto candidates_annotated.csv.

USAGE
-----
    python merge_psortb.py \
        --psortb PSORTb_results.txt \
        --candidates results/candidates_annotated.csv \
        --fasta webtools/all_clean.faa \
        --outdir results/

OUTPUT
------
  results/candidates_final.csv      full table, all filters joined
  results/surface_candidates.faa    OuterMembrane + Extracellular only
                                    <- the set that goes to Step 5/6
  results/localization_summary.csv  per-genome breakdown by localization

COLUMNS ADDED
-------------
  psortb_localization    final call (Extracellular / Periplasmic / Unknown / ...)
  psortb_score           0-10, PSORTb's own confidence
  psortb_modules         which analysis modules supported the call
  psortb_n_modules       how many did  <- read this before trusting the call
  score_Extracellular / score_OuterMembrane / score_Periplasmic /
  score_CytoplasmicMembrane / score_Cytoplasmic
  surface_exposed        yes / no

WHY psortb_n_modules MATTERS
----------------------------
PSORTb combines ~12 independent modules (SVMs, motif scans, BLAST against
proteins of known localization, signal peptide detection) and votes. A call
backed by one module is not the same as a call backed by four, even when the
final score looks identical. If every positive in your set traces back to the
SAME single module, that is a systematic pattern in that module -- not
independent evidence -- and it must be reported that way.

This matters especially for multiheme cytochromes: their amino acid
composition is extreme (very high Cys and His), which is exactly the kind of
signal a composition-based SVM can latch onto for the wrong reason.
"""

import argparse
import csv
import os
import re
import sys

LOCS = ["Extracellular", "OuterMembrane", "Periplasmic",
        "CytoplasmicMembrane", "Cytoplasmic"]
SURFACE = {"Extracellular", "OuterMembrane"}


def parse_psortb(path):
    """Parse PSORTb long-format output into {seq_id: {...}}."""
    text = open(path, encoding="utf-8", errors="replace").read().replace("\r", "")
    blocks = [b for b in re.split(r"^-{10,}$", text, flags=re.M) if "SeqID" in b]
    if not blocks:
        blocks = [b for b in text.split("---") if "SeqID" in b]
    if not blocks:
        sys.exit(f"ERROR: no PSORTb records found in {path}")

    out = {}
    for b in blocks:
        m_id = re.search(r"SeqID:\s*(\S+)", b)
        if not m_id:
            continue
        sid = m_id.group(1)

        # final prediction: either "Loc  score" or a bare "Unknown"
        m_fp = re.search(r"Final Prediction:\s*\n\s*(\S+)[ \t]*([\d.]*)", b)
        final = m_fp.group(1) if m_fp else ""
        score = m_fp.group(2) if m_fp and m_fp.group(2) else ""

        # all five localization scores
        scores = {}
        m_sec = re.search(r"Localization Scores:(.*?)Final Prediction:", b, re.S)
        if m_sec:
            for line in m_sec.group(1).strip().split("\n"):
                parts = line.split()
                if len(parts) == 2 and parts[0] in LOCS:
                    scores[parts[0]] = parts[1]

        # which analysis modules actually returned something
        mods = []
        m_rep = re.search(r"Analysis Report:(.*?)Localization Scores:", b, re.S)
        if m_rep:
            for line in m_rep.group(1).strip().split("\n"):
                parts = line.split()
                if len(parts) >= 2 and parts[1] != "Unknown":
                    mods.append(parts[0].rstrip("-"))

        row = {
            "psortb_localization": final,
            "psortb_score": score,
            "psortb_modules": ";".join(mods),
            "psortb_n_modules": str(len(mods)),
            "surface_exposed": "yes" if final in SURFACE else "no",
        }
        for loc in LOCS:
            row[f"score_{loc}"] = scores.get(loc, "")
        out[sid] = row
    return out


def read_fasta(path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--psortb", required=True,
                    help="PSORTb long-format output (repeat for several batches)",
                    nargs="+")
    ap.add_argument("--candidates", default="results/candidates_annotated.csv")
    ap.add_argument("--fasta", default="webtools/all_clean.faa")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    ps = {}
    for p in args.psortb:
        if not os.path.exists(p):
            sys.exit(f"ERROR: {p} not found")
        part = parse_psortb(p)
        print(f"{p}: {len(part)} records")
        ps.update(part)

    if not os.path.exists(args.candidates):
        sys.exit(f"ERROR: {args.candidates} not found")
    with open(args.candidates, newline="") as fh:
        cands = list(csv.DictReader(fh))

    print(f"PSORTb records : {len(ps)}")
    print(f"candidate rows : {len(cands)}")

    new_cols = (["psortb_localization", "psortb_score", "psortb_modules",
                 "psortb_n_modules", "surface_exposed"] +
                [f"score_{l}" for l in LOCS])

    cols = list(cands[0].keys())
    for c in new_cols:
        if c not in cols:
            cols.append(c)

    blank = {c: "" for c in new_cols}
    blank["psortb_localization"] = "not_submitted"
    blank["surface_exposed"] = "no"

    for row in cands:
        sid = row.get("seq_id", "")
        row.update(ps.get(sid, dict(blank)))

    out_csv = os.path.join(args.outdir, "candidates_final.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(cands)

    # --- surface-exposed FASTA ---
    keep = {r["seq_id"] for r in cands if r.get("surface_exposed") == "yes"}
    out_faa = os.path.join(args.outdir, "surface_candidates.faa")
    n = 0
    if os.path.exists(args.fasta):
        with open(out_faa, "w") as fh:
            for sid, seq in read_fasta(args.fasta):
                if sid in keep:
                    fh.write(f">{sid}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
                    n += 1

    # --- per-genome localization summary ---
    per = {}
    for r in cands:
        g = r.get("genome_id", "")
        d = per.setdefault(g, {
            "genome_id": g,
            "short_name": r.get("short_name", ""),
            "group": r.get("group", ""),
            "n_candidates": 0, "n_signalp_positive": 0, "n_surface": 0,
        })
        d["n_candidates"] += 1
        if r.get("has_signal_peptide") == "yes":
            d["n_signalp_positive"] += 1
        loc = r.get("psortb_localization", "")
        if loc and loc != "not_submitted":
            d[loc] = d.get(loc, 0) + 1
        if r.get("surface_exposed") == "yes":
            d["n_surface"] += 1

    seen_locs = sorted({k for d in per.values() for k in d
                        if k in LOCS or k == "Unknown"})
    fcols = (["genome_id", "short_name", "group", "n_candidates",
              "n_signalp_positive"] + seen_locs + ["n_surface"])
    out_sum = os.path.join(args.outdir, "localization_summary.csv")
    with open(out_sum, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fcols, extrasaction="ignore")
        w.writeheader()
        for d in sorted(per.values(), key=lambda x: -x["n_surface"]):
            w.writerow({c: d.get(c, 0) for c in fcols})

    # --- report ---
    from collections import Counter
    locs = Counter(r.get("psortb_localization", "?") for r in cands)
    print("\n" + "=" * 58)
    print("PSORTb localization")
    for k, v in locs.most_common():
        print(f"  {k:<22} {v:>4}")

    surf = [r for r in cands if r.get("surface_exposed") == "yes"]
    print(f"\nSurface-exposed (OM + Extracellular): {len(surf)}")
    if surf:
        modc = Counter(r.get("psortb_modules", "") for r in surf)
        print("  supporting-module patterns:")
        for k, v in modc.most_common():
            print(f"    {v:>3}x  {k or '(none)'}")
        single = sum(1 for r in surf if r.get("psortb_n_modules") == "1")
        if single:
            print(f"\n  ⚠️  {single}/{len(surf)} rest on a SINGLE module.")
        if len(modc) <= 2:
            print("  ⚠️  Nearly all positives share the same module pattern —"
                  "\n      report this as one line of evidence, not many.")

    print(f"\nwrote {out_csv}")
    print(f"wrote {out_faa}  ({n} sequences)")
    print(f"wrote {out_sum}")


if __name__ == "__main__":
    main()
