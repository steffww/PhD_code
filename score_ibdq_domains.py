#!/usr/bin/env python3
"""Append complete-case IBDQ-32 domain scores to the one-stop analysis bundle.

Standard item map:
 bowel 1,5,9,13,17,20,22,24,26,29
 systemic 2,6,10,14,18
 emotional 3,7,11,15,19,21,23,25,27,30,31,32
 social 4,8,12,16,28

Existing locked IBDQ_total is never overwritten. A domain score is missing
unless every item in that domain is an integer from 1 to 7. No imputation.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

BUNDLE = Path(__file__).resolve().parent
QUESTIONNAIRE = Path(os.environ.get(
    "ONE_IBD_QUESTIONNAIRE",
    "/Users/steven/claude_project/new_metagenomic_project/data/整合问卷20260831.xlsx",
))
FOLLOWUP = Path(os.environ.get(
    "ONE_IBD_FOLLOWUP",
    "/Users/steven/claude_project/new_metagenomic_project/data/Follow_20250831_xiangya.xlsx",
))
DOMAINS = {
    "IBDQ_bowel": [1, 5, 9, 13, 17, 20, 22, 24, 26, 29],
    "IBDQ_systemic": [2, 6, 10, 14, 18],
    "IBDQ_emotional": [3, 7, 11, 15, 19, 21, 23, 25, 27, 30, 31, 32],
    "IBDQ_social": [4, 8, 12, 16, 28],
}
assert sorted(x for values in DOMAINS.values() for x in values) == list(range(1, 33))
NEW = list(DOMAINS) + [
    "IBDQ_items_valid_n",
    "IBDQ_item_total_recomputed",
    "IBDQ_total_locked_discordant",
    "IBDQ_domain_scoring_status",
]

def score_items(frame: pd.DataFrame) -> pd.DataFrame:
    """frame has exactly 32 item columns in order; invalid values act missing."""
    assert frame.shape[1] == 32
    x = frame.apply(pd.to_numeric, errors="coerce")
    valid = x.notna() & x.ge(1) & x.le(7) & x.eq(np.floor(x))
    ans = pd.DataFrame(index=frame.index)
    for domain, items in DOMAINS.items():
        idx = [i - 1 for i in items]
        part = x.iloc[:, idx]
        ok = valid.iloc[:, idx].all(axis=1)
        ans[domain] = part.sum(axis=1, min_count=len(idx)).where(ok)
    ans["IBDQ_items_valid_n"] = valid.sum(axis=1)
    ans["IBDQ_item_total_recomputed"] = x.sum(axis=1, min_count=32).where(
        valid.all(axis=1)
    )
    ans["IBDQ_domain_scoring_status"] = np.where(
        valid.all(axis=1), "complete_32", "incomplete_or_invalid"
    )
    return ans

def main() -> None:
    for p in (QUESTIONNAIRE, FOLLOWUP):
        assert p.is_file(), p
    full_path = BUNDLE / "metadata" / "baseline_full_450.csv"
    core_path = BUNDLE / "metadata" / "baseline_core_450.csv"
    follow_path = BUNDLE / "metadata" / "followup_long_323.csv"
    full = pd.read_csv(full_path, dtype={"participant_id": str})
    core = pd.read_csv(core_path, dtype={"participant_id": str})
    visits = pd.read_csv(follow_path, dtype={"participant_id": str})
    assert len(full) == len(core) == 450 and len(visits) == 323
    assert set(full.sample_id) == set(core.sample_id)
    assert not set(NEW) & set(full.columns), "Domain columns already present; do not double-run."

    q = pd.read_excel(QUESTIONNAIRE, sheet_name="IBD问卷", dtype=str)
    q["样本编号"] = q["样本编号"].astype("string").str.strip()
    ibd = full[full.IBD.eq("IBD")]
    q = q[q["样本编号"].isin(set(ibd.participant_id))].copy()
    assert len(q) == 283 and q["样本编号"].is_unique
    # The actual questionnaire block is numbered 1..32 in this exact order.
    headers = list(q.columns[554:586])
    assert len(headers) == 32
    assert headers[0].startswith("1、") and headers[-1].startswith("32、")
    source = q.set_index("样本编号")
    baseline_scores = score_items(source[headers])
    baseline_scores.index.name = "participant_id"
    sample_map = ibd.set_index("participant_id")["sample_id"]
    assert set(source.index) == set(sample_map.index)
    assert (source.loc[sample_map.index, "宏基因组编号"].astype(str).str.strip().to_numpy()
            == sample_map.to_numpy()).all()
    source_total = pd.to_numeric(source["IBDQ总分"], errors="coerce")
    comparable = baseline_scores["IBDQ_item_total_recomputed"].notna() & source_total.notna()
    source_mismatch = (
        baseline_scores.loc[comparable, "IBDQ_item_total_recomputed"]
        != source_total.loc[comparable]
    )
    assert comparable.sum() == 277 and source_mismatch.sum() == 2

    for table, path in ((full, full_path), (core, core_path)):
        s = table["participant_id"].map(baseline_scores["IBDQ_item_total_recomputed"])
        for col in DOMAINS:
            table[col] = table["participant_id"].map(baseline_scores[col])
        table["IBDQ_items_valid_n"] = table["participant_id"].map(
            baseline_scores["IBDQ_items_valid_n"]
        )
        table["IBDQ_item_total_recomputed"] = s
        locked = pd.to_numeric(table["IBDQ_total"], errors="coerce")
        table["IBDQ_total_locked_discordant"] = np.where(
            s.notna() & locked.notna(), s.ne(locked), pd.NA
        )
        table["IBDQ_domain_scoring_status"] = table["participant_id"].map(
            baseline_scores["IBDQ_domain_scoring_status"]
        ).fillna("not_applicable_nonIBD")
        assert len(table) == 450 and table.sample_id.is_unique
        table.to_csv(path, index=False, na_rep="")

    f = pd.read_excel(FOLLOWUP, sheet_name="Sheet1", dtype=str)
    f["ID"] = f["ID"].astype("string").str.strip()
    f["source_excel_row"] = f.index + 2
    f = f[f["ID"].isin(set(visits.participant_id))].copy()
    fkey = f.set_index(["ID", "source_excel_row"])
    vkey = pd.MultiIndex.from_frame(
        visits[["participant_id", "source_excel_row"]].rename(
            columns={"participant_id": "ID"}
        )
    )
    assert vkey.is_unique and set(vkey).issubset(set(fkey.index))
    item_cols = [f"OutIBDQ{i}" for i in range(1, 33)]
    assert all(col in f.columns for col in item_cols)
    fscores = score_items(fkey[item_cols])
    fscores = fscores.loc[vkey].reset_index(drop=True)
    for col in DOMAINS:
        visits[col] = fscores[col].values
    visits["IBDQ_items_valid_n"] = fscores["IBDQ_items_valid_n"].values
    visits["IBDQ_item_total_recomputed"] = fscores["IBDQ_item_total_recomputed"].values
    prior = pd.to_numeric(visits["ibdq_followup"], errors="coerce")
    calc = visits["IBDQ_item_total_recomputed"]
    visits["IBDQ_total_locked_discordant"] = np.where(
        prior.notna() & calc.notna(), prior.ne(calc), pd.NA
    )
    visits["IBDQ_domain_scoring_status"] = fscores["IBDQ_domain_scoring_status"].values
    assert len(visits) == 323 and int(visits["IBDQ_total_locked_discordant"].eq(True).sum()) == 0
    visits.to_csv(follow_path, index=False, na_rep="")

    bfull = pd.read_csv(full_path)
    assert bfull.loc[bfull.IBD.eq("IBD"), "IBDQ_domain_scoring_status"].eq(
        "complete_32"
    ).sum() == 277
    print(
        "IBDQ_DOMAIN_SCORE_PASS",
        "baseline_complete=277/283",
        "locked_total_discordant=2",
        f"followup_complete={int(visits.IBDQ_item_total_recomputed.notna().sum())}/323",
    )

if __name__ == "__main__":
    main()
