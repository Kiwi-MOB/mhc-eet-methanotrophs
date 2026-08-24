# Extracellular electron transfer potential in aerobic methanotrophs and anammox bacteria

*An in silico comparative genomics survey*

> **Status:** in progress. Started July 2026.
> **Author:** CAI Wenjia ( Kiwi ), Dept. of Civil and Environmental Engineering, PolyU

---

## 1. Question

Coupled MOB–anammox systems are usually explained by **metabolite-mediated** interactions:
MOB supply organic carbon (methanol, acetate, formate) or reduce nitrate to nitrite, which
anammox then uses. A third possibility — **direct extracellular electron transfer (EET)**
between the two partners — is rarely examined, in part because nobody has asked whether
methanotrophs carry the genetic machinery for it.

**This project asks:** do aerobic methanotroph genomes encode candidate outer-membrane
multiheme cytochromes (MHCs) of the kind that mediate EET in known electroactive bacteria —
and do anammox genomes encode plausible counterparts on the receiving side?

This is a **potential** survey. It cannot demonstrate that EET occurs. It defines which
lineages and which proteins would be worth testing experimentally.

## 2. Approach

30 representative genomes spanning:

| Group | n | Rationale |
|---|---|---|
| Type I MOB (*Methylococcaceae*) | 10 | Includes *M. capsulatus* Bath, the strain used in our reactor work |
| Type II MOB (Alphaproteobacteria) | 4 | Phylogenetic breadth |
| Verrucomicrobial MOB | 2 | Deep outgroup among aerobic methanotrophs |
| Anaerobic methane oxidisers | 2 | *Ca.* Methylomirabilis, *Ca.* Methanoperedens — the latter with documented EET machinery |
| Anammox (*Brocadiaceae*) | 6 | The partner side of the couple |
| Controls | 6 | *Geobacter*, *Shewanella* (positive); *E. coli*, *P. aeruginosa* (negative); *Nitrosomonas*, *Nitrospira* (methodological / ecological) |

Full accession list and selection rationale: [`genomes.csv`](genomes.csv)

### Pipeline

```
proteomes/*.faa
      |
      |  scan_mhc.py          CXXCH / CXXCK motif detection (overlap-aware)
      v
candidates.faa  +  candidates.csv        [>= 3 heme motifs]
      |
      |  SignalP 6.0          signal peptide
      |  DeepTMHMM            transmembrane topology
      |  PSORTb               subcellular localisation
      |  InterProScan         domain architecture
      v
filtered candidates                       [secreted / outer membrane]
      |
      |  NCBI gene neighbourhood inspection    porin-cytochrome operon?
      |  AlphaFold3 (with HEM ligands)         fold and heme packing
      v
final candidate set + comparative figures
```

### Scoring criteria (defined before analysis)

A protein is a **high-confidence EET candidate** if it meets all of:

1. ≥ 3 CXXCH/CXXCK heme-binding motifs
2. Predicted signal peptide (Sec or Tat)
3. Predicted outer-membrane or extracellular localisation
4. Annotated cytochrome c domain (InterPro/Pfam)

**Supporting evidence** (not required): a porin-family gene within 5 genes of the candidate.

*These criteria were fixed in advance and are not adjusted post hoc.*

## 3. Repository contents

```
.
├── README.md
├── genomes.csv                 genome selection + accessions + rationale
├── scan_mhc.py                 motif scanner (no dependencies)
├── proteomes/                  input .faa files (not tracked; see below)
├── results/
│   ├── per_genome_summary.csv  MHC inventory per genome
│   ├── candidates.csv          master candidate table
│   └── candidates.faa          candidate sequences
└── figures/
```

Proteome FASTA files are not committed. Reproduce them by downloading the accessions
in `genomes.csv` from NCBI, naming each file `<accession>.faa`.

## 4. Reproducing

```bash
python scan_mhc.py --indir proteomes/ --outdir results/ --min-hemes 3
```

Requires Python 3.8+. No external packages.

## 5. Findings so far

### Localisation filtering (SignalP 6.0 → PSORTb 3.0)

Of 374 multiheme candidates, 239 carried a predicted signal peptide, and 
PSORTb assigned 15 to an extracellular localisation. No candidate in the 
entire dataset was assigned to the outer membrane.

**Positive control behaviour was asymmetric, and this constrains 
interpretation.** PSORTb correctly recovered OmcS, an experimentally 
confirmed extracellular cytochrome of G. sulfurreducens. However, it 
assigned all S. oneidensis candidates to the periplasm or to Unknown, 
including the MtrC/OmcA system — outer-membrane lipoproteins that are 
experimentally confirmed to be surface-exposed. These are neither freely 
secreted nor beta-barrel proteins, and fall outside the structural classes 
PSORTb resolves.

Consequently, in this dataset a PSORTb "Extracellular" call carries 
information, but a "Periplasmic" or "Unknown" call does NOT exclude surface 
exposure. M. capsulatus Bath returned no surface-exposed candidates, but its 
profile is indistinguishable from that of S. oneidensis, where the same 
result is a known false negative. Bath's surface localisation therefore 
remains undetermined rather than negative, and is being addressed through 
gene-neighbourhood analysis and lipoprotein sorting rules.

All 15 extracellular calls were driven by the ECSVM module, with independent 
BLAST support in only one case. Multiheme cytochromes have atypical amino 
acid composition (high Cys and His), which a composition-based classifier may 
respond to for reasons unrelated to localisation. These 15 calls should be 
treated as one line of evidence, not fifteen.

## 6. Known limitations

- CXXCH motif counting detects the canonical heme attachment site only; rarer variants
  (CXXXCH, CX15CH) are missed.
- Localisation prediction for Planctomycetes (anammox) is unreliable — their cell plan
  differs from the Gram-negative model these predictors were trained on. Anammox results
  are therefore reported with lower confidence than the MOB results.
- Genetic potential is not activity. Nothing here shows these proteins are expressed,
  correctly matured, or functional.
- Genome selection is representative, not exhaustive; abundance in real reactors is not
  addressed.

## 7. License

Code: MIT. Analysis outputs: CC-BY-4.0.
