#!/usr/bin/env python3
"""
erepo_combination_benchmark.py

Re-classifies every ClinGen Evidence Repository (ERepo) interpretation from the
expert panel's own applied ("Met") evidence codes, using three published combiners:

  R15   ACMG/AMP 2015 combining rules (Richards et al. 2015, Table 5)
  R15S  R15 + very strong pathogenic evidence with >= 1 supporting pathogenic
        criterion = Likely pathogenic (generalizes ClinGen SVI PVS1 + PM2_Supporting)
  PTS   Naturally scaled point system (Tavtigian et al. 2018, 2020):
        supporting 1, moderate 2, strong 4, very strong / stand-alone 8;
        P >= 10, LP 6..9, VUS 0..5, LB -1..-6, B <= -7

Usage:
  curl -o erepo.tsv "https://erepo.clinicalgenome.org/evrepo/api/classifications/all?format=tabbed"
  python3 erepo_combination_benchmark.py erepo.tsv out.tsv

Requires: Python >= 3.9, pandas. Licence: MIT.
"""
import re
import sys
import pandas as pd

DEFAULT = {"PVS": "very_strong", "PS": "strong", "PM": "moderate", "PP": "supporting",
           "BA": "stand_alone", "BS": "strong", "BP": "supporting"}
SUFFIX = {"supporting": "supporting", "moderate": "moderate", "strong": "strong",
          "very strong": "very_strong", "stand alone": "stand_alone"}
POINTS = {"supporting": 1, "moderate": 2, "strong": 4, "very_strong": 8, "stand_alone": 8}
ASSERTION = {"pathogenic": "P", "likely pathogenic": "LP", "uncertain significance": "VUS",
             "likely benign": "LB", "benign": "B"}


def parse_codes(field):
    """'PVS1, PM2_Supporting, PP4_Moderate' -> [(code, strength, direction), ...]"""
    out = []
    for tok in (t.strip() for t in field.split(",")):
        if not tok:
            continue
        if "_" in tok:
            code, suf = tok.split("_", 1)
            strength = SUFFIX[suf.strip().lower()]
        else:
            code = tok
            strength = DEFAULT[re.match(r"(PVS|PS|PM|PP|BA|BS|BP)\d", code).group(1)]
        out.append((code.upper(), strength, "benign" if code[0] == "B" else "pathogenic"))
    return out


def rules_2015(codes, svi=False):
    vs = s = m = p = sa = bs = bp = 0
    for _, st, d in codes:
        if d == "pathogenic":
            if st == "very_strong":
                vs += 1
            elif st == "strong":
                s += 1
            elif st == "moderate":
                m += 1
            else:
                p += 1
        else:
            if st == "stand_alone":
                sa += 1
            elif st in ("strong", "very_strong"):
                bs += 1
            else:
                bp += 1
    if vs > 1:                      # additional very-strong evidence counted as strong
        s, vs = s + vs - 1, 1
    is_p = ((vs and (s >= 1 or m >= 2 or (m == 1 and p >= 1) or p >= 2)) or s >= 2
            or (s == 1 and (m >= 3 or (m == 2 and p >= 2) or (m == 1 and p >= 4))))
    is_lp = ((vs and m >= 1) or (s == 1 and 1 <= m <= 2) or (s == 1 and p >= 2) or m >= 3
             or (m == 2 and p >= 2) or (m == 1 and p >= 4) or (svi and vs and p >= 1))
    is_b = sa >= 1 or bs >= 2
    is_lb = (bs == 1 and bp >= 1) or bp >= 2
    path = "P" if is_p else ("LP" if is_lp else None)
    ben = "B" if is_b else ("LB" if is_lb else None)
    if path and ben:                # contradictory evidence
        return "VUS"
    return path or ben or "VUS"


def point_sum(codes):
    return sum(POINTS[st] * (1 if d == "pathogenic" else -1) for _, st, d in codes)


def points_class(score):
    if score >= 10:
        return "P"
    if score >= 6:
        return "LP"
    if score >= 0:
        return "VUS"
    if score >= -6:
        return "LB"
    return "B"


def main(snapshot, out_tsv):
    df = pd.read_csv(snapshot, sep="\t", dtype=str, keep_default_na=False)
    df = df[df["Retracted"] != "true"]
    rows = []
    for _, r in df.iterrows():
        codes = parse_codes(r["Applied Evidence Codes (Met)"])
        if not codes:
            continue
        rows.append({"uuid": r["Uuid"], "vcep": " ".join(r["Expert Panel"].split()),
                     "gene": r["HGNC Gene Symbol"], "met_codes": r["Applied Evidence Codes (Met)"],
                     "assertion": ASSERTION[r["Assertion"].strip().lower()],
                     "points": point_sum(codes), "R15": rules_2015(codes),
                     "R15S": rules_2015(codes, svi=True), "PTS": points_class(point_sum(codes))})
    res = pd.DataFrame(rows)
    res.to_csv(out_tsv, sep="\t", index=False)
    n = len(res)
    for m in ["R15", "R15S", "PTS"]:
        k = int((res[m] == res.assertion).sum())
        print(f"{m:5s} exact agreement {k}/{n} = {100 * k / n:.2f}%")


if __name__ == "__main__":
    main(*sys.argv[1:3])
