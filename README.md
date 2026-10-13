# ERepo combination benchmark

Code to reproduce the published-combiner analyses of *Points versus rules: reproducing 13,012 ClinGen expert-panel variant classifications from their own evidence codes* (Jiménez Nieto, 2026).

Every interpretation in the ClinGen Evidence Repository (ERepo) is re-classified from the expert panel's own applied ("Met") evidence codes and strengths with three published combiners, and compared with the panel's assertion:

| Combiner | Logic |
|---|---|
| `R15` | ACMG/AMP 2015 combining rules (Richards et al. 2015, Table 5) |
| `R15S` | `R15` + very strong pathogenic evidence (PVS1 or an upgraded criterion) with ≥1 supporting pathogenic criterion → Likely pathogenic (generalizing the ClinGen SVI PM2 v1.0 recommendation) |
| `PTS` | Naturally scaled point system (Tavtigian et al. 2020): 1/2/4/8 points; P ≥ 10, LP 6–9, VUS 0–5, LB −1 to −6, B ≤ −7 |

## Reproduce

```bash
curl -o erepo.tsv "https://erepo.clinicalgenome.org/evrepo/api/classifications/all?format=tabbed"
python3 erepo_combination_benchmark.py erepo.tsv results.tsv
```

The study snapshot was downloaded on 6 October 2026 (SHA-256 `2adc311c013928a8ebb0a94e749cc2557d868ff89771d300551a6dbd32d83dbf`, 13,293 rows). ERepo is continually updated, so a fresh download will contain more interpretations; on the study snapshot the script prints:

```
R15   exact agreement 12012/13012 = 92.31%
R15S  exact agreement 12377/13012 = 95.12%
PTS   exact agreement 12348/13012 = 94.90%
```

## Files

- `erepo_combination_benchmark.py` — the three combiners and the evaluation (Python ≥ 3.9, pandas).

## Licence

Code: MIT. ERepo data are provided by ClinGen; please cite ClinGen and the relevant expert panels when reusing them.

## Citation

If you use this code, please cite the article (Jiménez Nieto N, 2026; full reference to be added on publication) and the archived software release (see `CITATION.cff`).
