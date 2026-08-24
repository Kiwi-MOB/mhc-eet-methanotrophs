#!/usr/bin/env python3
"""
merge_signalp.py
----------------
Join SignalP 6.0 predictions back onto candidates.csv, and write out the
signal-peptide-positive subset for the next annotation round.

USAGE
-----
    python merge_signalp.py \
        --signalp prediction_results.txt \
        --idmap webtools/id_map.csv \
        --candidates results/candidates.csv \
        --fasta webtools/all_clean.faa \
        --metadata genomes.csv \
        --outdir results/

OUTPUT
------
  results/candidates_annotated.csv   candidates.csv + SignalP columns
  results/signalp_positive.faa       only sequences WITH a signal peptide
                                     <- submit this to PSORTb / DeepTMHMM
  results/funnel_summary.csv         per-genome candidates -> SP-positive

NEW COLUMNS ADDED
-----------------
  seq_id              the seqNNNN handle used with the web tools
  signalp_prediction  SP / LIPO / TAT / TATLIPO / PILIN / OTHER
  signalp_prob        probability of the winning class
  signalp_confidence  high (>=0.90) / medium (0.70-0.90) / low (<0.70)
  has_signal_peptide  yes / no
  cs_position         cleavage site, as reported

WHY signalp_confidence EXISTS
-----------------------------
SignalP reports a hard call, but that call is just whichever class scored
highest. A sequence can be labelled "SP" on a score of 0.51 against 0.48 for
OTHER -- effectively a coin flip. Those calls should not be treated the same
as a 0.999 call. Anything below 0.70 is flagged so it can be handled
separately (or excluded) rather than silently counted as a positive.
"""

import argparse
import csv
import os
import sys

CONF_HIGH = 0.90
CONF_MED = 0.70


def read_signalp(path):
    """
    Parse a SignalP 6.0 'Prediction summary' file.

    Format (tab-separated, '#' comment lines):
      # SignalP-6.0  Organism: Other  Timestamp: ...
      # ID  Prediction  OTHER  SP(Sec/SPI)  LIPO(Sec/SPII)  ...  CS Position
      seq0001  SP  0.484538  0.512940  ...  CS pos: 40-41. Pr: 0.4893
    """
    header = None
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("#"):
                fields = [c.strip() for c in line.lstrip("#").split("\t")]
                # the column-name line is the one starting with "ID"
                if fields and fields[0].upper() == "ID":
                    header = fields
                continue
            if header is None:
                sys.exit("ERROR: no '# ID ...' column header found in "
                         f"{path}. Is this the 'Prediction summary' file?")

            parts = line.split("\t")
            row = dict(zip(header, parts))          # ragged rows are fine
            sid = parts[0]

            # find the winning probability column, whatever it is called
            best_col, best_val = None, -1.0
            for col in header[2:]:
                if col.lower().startswith("cs"):
                    continue
                try:
                    v = float(row.get(col, ""))
                except (TypeError, ValueError):
                    continue
                if v > best_val:
                    best_col, best_val = col, v

            has_sp = bool(best_col) and not best_col.upper().startswith("OTHER")
            if best_val >= CONF_HIGH:
                conf = "high"
            elif best_val >= CONF_MED:
                conf = "medium"
            else:
                conf = "low"

            out[sid] = {
                "signalp_prediction": row.get("Prediction", ""),
                "signalp_prob": f"{best_val:.4f}" if best_val >= 0 else "",
                "signalp_confidence": conf,
                "has_signal_peptide": "yes" if has_sp else "no",
                "cs_position": row.get("CS Position", "").strip(),
            }
    return out


def read_csv_rows(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


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
    ap.add_argument("--signalp", required=True)
    ap.add_argument("--idmap", default="webtools/id_map.csv")
    ap.add_argument("--candidates", default="results/candidates.csv")
    ap.add_argument("--fasta", default="webtools/all_clean.faa")
    ap.add_argument("--metadata", default="genomes.csv")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    for p in (args.signalp, args.idmap, args.candidates):
        if not os.path.exists(p):
            sys.exit(f"ERROR: {p} not found")

    os.makedirs(args.outdir, exist_ok=True)

    sp = read_signalp(args.signalp)
    idmap = read_csv_rows(args.idmap)
    cands = read_csv_rows(args.candidates)
    meta = {r["genome_id"]: r for r in read_csv_rows(args.metadata)
            if r.get("genome_id")}

    print(f"SignalP rows   : {len(sp)}")
    print(f"id_map rows    : {len(idmap)}")
    print(f"candidate rows : {len(cands)}")
    if meta:
        print(f"genomes.csv    : {len(meta)} genomes")
    else:
        print("genomes.csv    : NOT FOUND -- organism names will be missing")

    # seq_id lookup keyed on (genome_id, protein_acc)
    key2seq = {(r["genome_id"], r["protein_acc"]): r["seq_id"] for r in idmap}

    meta_cols = []
    if meta:
        sample = next(iter(meta.values()))
        meta_cols = [c for c in ("organism", "short_name", "group")
                     if c in sample]

    new_cols = ["seq_id", "signalp_prediction", "signalp_prob",
                "signalp_confidence", "has_signal_peptide", "cs_position"]

    base_cols = [c for c in cands[0].keys()] if cands else []
    for c in meta_cols + new_cols:
        if c not in base_cols:
            base_cols.append(c)

    merged, unmatched = [], 0
    for row in cands:
        key = (row.get("genome_id", ""), row.get("protein_acc", ""))
        sid = key2seq.get(key, "")
        row["seq_id"] = sid
        if sid and sid in sp:
            row.update(sp[sid])
        else:
            unmatched += 1
            row.update({c: "" for c in new_cols[1:]})
        for c in meta_cols:
            if not row.get(c):
                row[c] = meta.get(row.get("genome_id", ""), {}).get(c, "")
        merged.append(row)

    if unmatched:
        print(f"\n⚠️  {unmatched} candidates had no SignalP result "
              f"(check that all sequences were submitted)")

    out_csv = os.path.join(args.outdir, "candidates_annotated.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=base_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(merged)

    # --- positive subset FASTA ---
    keep = {r["seq_id"] for r in merged if r.get("has_signal_peptide") == "yes"}
    out_faa = os.path.join(args.outdir, "signalp_positive.faa")
    n_written = 0
    if os.path.exists(args.fasta):
        with open(out_faa, "w") as fh:
            for sid, seq in read_fasta(args.fasta):
                if sid in keep:
                    fh.write(f">{sid}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
                    n_written += 1

    # --- funnel table ---
    funnel = {}
    for r in merged:
        g = r.get("genome_id", "")
        d = funnel.setdefault(g, {
            "genome_id": g,
            "short_name": r.get("short_name", ""),
            "group": r.get("group", ""),
            "n_candidates": 0, "n_signalp_positive": 0,
            "n_high_confidence": 0,
        })
        d["n_candidates"] += 1
        if r.get("has_signal_peptide") == "yes":
            d["n_signalp_positive"] += 1
            if r.get("signalp_confidence") == "high":
                d["n_high_confidence"] += 1

    out_funnel = os.path.join(args.outdir, "funnel_summary.csv")
    fcols = ["genome_id", "short_name", "group", "n_candidates",
             "n_signalp_positive", "n_high_confidence"]
    with open(out_funnel, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fcols, extrasaction="ignore")
        w.writeheader()
        for d in sorted(funnel.values(),
                        key=lambda x: -x["n_signalp_positive"]):
            w.writerow(d)

    # --- report ---
    types, confs = {}, {}
    for r in merged:
        types[r.get("signalp_prediction", "?")] = types.get(
            r.get("signalp_prediction", "?"), 0) + 1
        if r.get("has_signal_peptide") == "yes":
            c = r.get("signalp_confidence", "?")
            confs[c] = confs.get(c, 0) + 1

    n_pos = sum(1 for r in merged if r.get("has_signal_peptide") == "yes")
    print("\n" + "=" * 58)
    print("SignalP prediction classes")
    for k, v in sorted(types.items(), key=lambda x: -x[1]):
        print(f"  {k:<12} {v:>5}")
    print(f"\nWith a signal peptide : {n_pos} / {len(merged)}"
          f"  ({100*n_pos/max(len(merged),1):.0f}%)")
    print("  confidence breakdown:")
    for k in ("high", "medium", "low"):
        if k in confs:
            print(f"    {k:<8} {confs[k]:>5}")
    if confs.get("low"):
        print(f"\n  ⚠️  {confs['low']} positives scored below {CONF_MED} —"
              " treat these as uncertain, not as findings.")

    print(f"\nwrote {out_csv}")
    print(f"wrote {out_faa}  ({n_written} sequences)")
    print(f"wrote {out_funnel}")
    print("\nNext: submit signalp_positive.faa to PSORTb and DeepTMHMM.")


if __name__ == "__main__":
    main()
