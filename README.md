# Genome-based prediction of extracellular electron transfer fails against a proteomic benchmark

*A case study in* Methylococcus capsulatus *(Bath)*

> **Status:** benchmark complete. **Author:** CAI Wenjia (Kiwi), Dept. of
> Civil and Environmental Engineering, The Hong Kong Polytechnic University.
> **Context:** independent computational project supporting laboratory research
> on anammox × methane-oxidising bacteria (MOB) coupling.

---

## 1. Where this started

Our laboratory observes coupling between methane-oxidising bacteria and anammox
bacteria, with dissimilatory nitrate reduction to ammonium (DNRA) occurring in
the coupled condition. Showing *that* coupling occurs does not explain *how*.
Three mechanisms are plausible:

| | Mechanism | What passes between partners |
|---|---|---|
| H1 | Metabolite-mediated | methanol, formate, acetate |
| H2 | Nitrogen-mediated | nitrite or ammonium |
| **H3** | **Direct interspecies electron transfer** | **electrons, cell surface to cell surface** |

H3 requires hardware: surface-exposed multiheme *c*-type cytochromes, typically
encoded alongside an outer-membrane porin forming a conduit. Whether aerobic
methanotrophs carry this hardware looked answerable from sequence alone, so I
built a four-stage genomic screen to find out.

**Then I found the screen could be tested.** *M. capsulatus* Bath has a
proteomic gold standard: Larsen & Karlsen (2016) identified, by LTQ-ORBITRAP
mass spectrometry of surface-extracted fractions, eight copper-responsive
*c*-type cytochromes physically located at the cell surface. Tanaka et al.
(2018) then demonstrated EET in Bath functionally, with knockout of one of
these genes (MCA0421) significantly suppressing anodic current.

So before using the screen to answer anything, I measured what it recovers.

## 2. Benchmark result

**Recall: 0/8.**

| Locus tag | Hemes | L1 motif ≥3 | L2 SignalP | L3 PSORTb | Surface called |
|---|---|---|---|---|---|
| MCA0421 | 9 | ✅ | SP | Periplasmic | ❌ |
| MCA0423 | 8 | ✅ | SP | Unknown | ❌ |
| MCA0424 | 7 | ✅ | SP | **Cytoplasmic** | ❌ |
| MCA0426 | 7 | ✅ | SP | Unknown | ❌ |
| MCA2259 | 9 | ✅ | SP | Periplasmic | ❌ |
| MCA0444 | <3 | ❌ | — | — | ❌ |
| MCA1906 | <3 | ❌ | — | — | ❌ |
| MCA2590 | <3 | ❌ | — | — | ❌ |

| Stage | Recall |
|---|---|
| L1 — heme motif count ≥3 | **5/8** |
| L2 — signal peptide (SignalP 6.0) | **5/5 of those reaching it** |
| L3 — localisation (PSORTb 3.0) | **0/5** |
| **End to end** | **0/8** |

The screen called **zero** Bath proteins surface-exposed. Precision is
undefined because there were no positive calls to evaluate.

## 3. Where it fails, and why

**The failure is not uniform.** Signal peptide prediction was lossless on this
set (5/5). Two distinct mechanisms account for the rest.

### Failure 1 — motif threshold excludes diheme cytochromes

MCA0444, MCA1906 and MCA2590 are MauG-family proteins, which are diheme. A
threshold of ≥3 heme motifs — chosen to suppress chance matches — removes an
entire structural class that is genuinely surface-located. **Three of eight
positives are lost before any localisation step runs.**

### Failure 2 — localisation prediction is not merely uncertain, it is wrong

PSORTb misassigned all five remaining positives, and not in a consistent
direction: two to the periplasm, two to Unknown, and **MCA0424 to the
cytoplasm** — a protein with a predicted signal peptide, seven heme motifs,
and direct mass-spectrometric detection at the cell surface.

This is consistent with what the same tool did to the positive control:
PSORTb assigned every member of the *S. oneidensis* mtr operon
(MtrA/MtrC/OmcA/MtrF, WP_011071900–903) to Unknown or Periplasmic, although
RefSeq itself annotates MtrC as an "extracellular … surface decaheme
cytochrome". Outer-membrane lipoproteins are neither freely secreted nor
beta-barrels and fall outside the structural classes PSORTb resolves.

### Failure 3 — the porin criterion fails on annotation, not algorithm

Manual gene-neighbourhood inspection of five Bath candidates found no adjacent
porin, which I initially read as evidence against a conduit. **That reading was
wrong.** Larsen & Karlsen report that MCA0427, a beta-barrel outer membrane
protein, sits within this same operon — an eleven-gene cluster containing eight
multiheme cytochromes with 58 CxxCH motifs in total — and note that in
*Shewanella* the equivalent beta-barrel MtrB forms the porin bridging
periplasmic MtrA and surface-exposed MtrC.

The porin is there. NCBI does not annotate it as one, so a criterion based on
reading annotation text cannot see it. **This failure lies in database
completeness, not in the prediction algorithm** — a second, independent route
to the same wrong conclusion.

## 4. Conclusions

**On the method.** On the only aerobic methanotroph with a proteomic gold
standard, sequence-based prediction of surface-exposed multiheme cytochromes
recovered none of eight validated positives. Failure occurred at two
independent stages, by two unrelated mechanisms, with a third failure arising
from annotation rather than prediction. Negative results from this class of
pipeline should not be interpreted as absence of EET capability.

**On my own scan.** This repository also contains a screen of 28 genomes
spanning aerobic and anaerobic methane oxidisers, anammox bacteria and
controls (`results/`). Given the benchmark above, **those results cannot
support conclusions about the EET capability of any taxon**, including the
question that motivated the project. They are retained as the substrate of the
benchmark and as a description of multiheme cytochrome inventories, not as
evidence of absence.

**On the coupling question.** DIET is not excluded — it is *untested* by
genomic means, because genomic means do not work here. The coupling observed
in our reactors must be discriminated experimentally.

## 5. Recommendations

For anyone screening non-model organisms for EET machinery:

1. **Set the heme motif threshold at ≥2, not ≥3.** Diheme MauG-family
   cytochromes are surface-located and are lost at ≥3. The chance-match cost of
   ≥2 is lower than the cost of removing a structural class.
2. **Do not use a negative PSORTb call to exclude surface exposure.** In this
   dataset it was wrong 5/5 on validated positives and 0/4 on the *Shewanella*
   mtr operon.
3. **Detect porins by domain search, not annotation text.** MCA0427 is a
   beta-barrel outer membrane protein that no keyword search on NCBI
   annotation will find. Use Pfam beta-barrel HMMs against the neighbourhood.
4. **Benchmark before concluding.** Both organisms in this study with
   experimentally validated EET — *S. oneidensis* and *M. capsulatus* Bath —
   returned false negatives. Any pipeline of this type should be run against a
   validated positive set before its output is interpreted.

## 6. Data and methods

```
proteomes (28 genomes, ~106,000 proteins)
   │  scan_mhc.py — CXXCH / CXXCK motifs, ≥3
   ↓  374 candidates
   │  SignalP 6.0
   ↓  239 with signal peptide
   │  PSORTb 3.0
   ↓  15 called extracellular (none in Bath)
   │  manual gene neighbourhood inspection
```

| File | Contents |
|---|---|
| `results/benchmark_recall.csv` | **the benchmark: 8 gold-standard proteins, stage by stage** |
| `genomes.csv` | genome selection, accessions, quality metrics |
| `results/candidates_final.csv` | all 374 candidates, all filters joined |
| `results/localization_summary.csv` | per-genome localisation breakdown |
| `results/gene_neighbourhood.csv` | manual locus inspection |

Scripts: `scan_mhc.py`, `prep_for_webtools.py`, `merge_signalp.py`,
`merge_psortb.py`. Python 3.8+, no external packages. Proteome FASTA files are
not committed — download the accessions in `genomes.csv` from NCBI.

## 7. Limitations of the benchmark itself

- **n = 8, one organism.** This is a case study, not a general performance
  estimate. Whether recall is this poor across taxa is unknown.
- The gold standard is mass-spectrometric detection in a surface-enriched
  fraction; enrichment is imperfect and some assignments may not be strictly
  outer-surface.
- Only one localisation predictor was tested. Others (BUSCA, CELLO,
  DeepLocPro) may perform differently and were not evaluated.
- Gene neighbourhoods were inspected manually for six proteins, not
  systematically.
- Bath's multiheme cytochromes are copper-repressed; the gold standard derives
  from copper-limited cells. Genomic prediction is condition-blind by nature,
  but this means the benchmark targets are proteins expressed under one
  specific regime.

## 8. References

- Larsen Ø & Karlsen OA (2016) *MicrobiologyOpen* 5:254–267. doi:10.1002/mbo3.324
- Tanaka K *et al.* (2018) *Front Microbiol* 9:2905. doi:10.3389/fmicb.2018.02905
- Ward N *et al.* (2004) *PLoS Biol* 2:e303.

## 9. Tooling note

Analysis scripts were written with AI assistance (Claude). Study design, genome
curation, benchmark construction, literature tracing, gene neighbourhood
inspection, and interpretation are my own.

## 10. License

Code: MIT. Analysis outputs: CC-BY-4.0.
