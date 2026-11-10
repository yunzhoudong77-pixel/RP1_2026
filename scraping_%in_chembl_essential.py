#!/usr/bin/env python3
"""Survey single-dose-point (% inhibition) coverage in a local ChEMBL SQLite release.

Replaces the earlier REST-API version. Querying the local database removes the
rate limiting and, more importantly, makes it possible to split document counts
by data source, which the API could not do in a single pass.

Usage
-----
    python survey_local.py --db /path/to/chembl_37.db
    python survey_local.py --db ... --year 2019
    python survey_local.py --db ... --schema      # inspect tables and exit

Outputs (written next to the script)
------------------------------------
    survey_by_type.csv        activity counts per target per standard_type
    survey_inhibition.csv     the main per-target % inhibition breakdown
    inhibition_documents.csv  one row per (target, document): candidate papers
                              to check supplementary information against
    inhibition_assays.csv     distinct assay descriptions, where the test
                              concentration and any panel name are recorded
"""

from __future__ import annotations

import argparse
import csv
import os
import sqlite3
import sys

# ChEMBL source id 1 is "Scientific Literature"; all other ids are datasets
# deposited in bulk (SGC, EUbOPEN/PDSP, PubChem, ...).
SRC_LITERATURE = 1

STANDARD_TYPES = ["IC50", "Ki", "Kd", "EC50", "AC50",
                  "Inhibition", "Activity", "Potency"]

# chembl_id, label, gene  — the 32 "essential" anti-targets of the avoidome,
# de-duplicated (the figure lists SLCO1B1 twice, once as "SLCOB1" and once as
# the protein name OATP1B1).
TARGETS: list[tuple[str, str, str]] = [
    ("CHEMBL4302",    "ABCB1",   "ABCB1"),
    ("CHEMBL5393",    "BCRP",    "ABCG2"),
    ("CHEMBL6020",    "ABCB11",  "ABCB11"),
    ("CHEMBL5685",    "SLC22A1", "SLC22A1"),
    ("CHEMBL1743122", "SLC22A2", "SLC22A2"),
    ("CHEMBL1641347", "SLC22A6", "SLC22A6"),
    ("CHEMBL3004",    "ABCC1",   "ABCC1"),
    ("CHEMBL1641348", "SLC22A8", "SLC22A8"),
    ("CHEMBL1697668", "SLCO1B1", "SLCO1B1"),
    ("CHEMBL1743121", "SLCO1B3", "SLCO1B3"),
    ("CHEMBL1743126", "SLC47A1", "SLC47A1"),
    ("CHEMBL1743127", "SLC47A2", "SLC47A2"),
    ("CHEMBL3253",    "HSA",     "ALB"),
    ("CHEMBL3356",    "CYP1A2",  "CYP1A2"),
    ("CHEMBL289",     "CYP2D6",  "CYP2D6"),
    ("CHEMBL3397",    "CYP2C9",  "CYP2C9"),
    ("CHEMBL340",     "CYP3A4",  "CYP3A4"),
    ("CHEMBL3401",    "PXR",     "NR1I2"),
    ("CHEMBL3201",    "AhR",     "AHR"),
    ("CHEMBL216",     "M1",      "CHRM1"),
    ("CHEMBL211",     "M2",      "CHRM2"),
    ("CHEMBL245",     "M3",      "CHRM3"),
    ("CHEMBL226",     "A1AR",    "ADORA1"),
    ("CHEMBL213",     "B1AR",    "ADRB1"),
    ("CHEMBL240",     "hERG",    "KCNH2"),
    ("CHEMBL1940",    "CaV1.2",  "CACNA1C"),
    ("CHEMBL231",     "H1",      "HRH1"),
    ("CHEMBL1941",    "H2",      "HRH2"),
    ("CHEMBL1833",    "5-HT2B",  "HTR2B"),
    ("CHEMBL251",     "A2AR",    "ADORA2A"),
    ("CHEMBL210",     "B2AR",    "ADRB2"),
    ("CHEMBL1980",    "NaV1.5",  "SCN5A"),
]


# --------------------------------------------------------------------------- #
# schema handling
# --------------------------------------------------------------------------- #
def columns(con: sqlite3.Connection, table: str) -> set[str]:
    return {r[1] for r in con.execute(f"PRAGMA table_info({table})")}


def locate_src(con: sqlite3.Connection) -> str:
    """Return the SQL expression giving a record's source id.

    `src_id` lives on `activities` in recent ChEMBL releases but has sat on
    `assays` and on `docs` in others, so resolve it rather than assume.
    """
    for table, alias in (("activities", "act"), ("assays", "a"), ("docs", "d")):
        if "src_id" in columns(con, table):
            print(f"  source id taken from {table}.src_id")
            return f"{alias}.src_id"
    sys.exit("No src_id column found on activities, assays or docs.")


def show_schema(con: sqlite3.Connection) -> None:
    for t in ("activities", "assays", "docs", "target_dictionary",
              "molecule_dictionary", "compound_structures", "source"):
        try:
            cols = sorted(columns(con, t))
        except sqlite3.Error:
            print(f"{t}: MISSING")
            continue
        print(f"\n{t} ({len(cols)} columns)")
        print("  " + ", ".join(cols))
    v = con.execute("SELECT name FROM sqlite_master WHERE type='table' "
                    "AND name LIKE 'version%'").fetchall()
    if v:
        for row in con.execute(f"SELECT * FROM {v[0][0]}"):
            print("\nversion:", row)


# --------------------------------------------------------------------------- #
# queries
# --------------------------------------------------------------------------- #
BASE_JOIN = """
    FROM activities act
    JOIN assays a             ON a.assay_id = act.assay_id
    JOIN target_dictionary td ON td.tid     = a.tid
    JOIN docs d               ON d.doc_id   = act.doc_id
"""


def by_type(con, ids, year):
    sql = f"""
        SELECT td.chembl_id                                        AS target,
               COALESCE(act.standard_type, '(none)')               AS std_type,
               COUNT(*)                                            AS n_all,
               SUM(CASE WHEN d.year >= ? THEN 1 ELSE 0 END)        AS n_recent
        {BASE_JOIN}
        WHERE td.chembl_id IN ({','.join('?' * len(ids))})
        GROUP BY td.chembl_id, std_type
    """
    return con.execute(sql, [year, *ids]).fetchall()


def inhibition_summary(con, ids, year, src):
    """Per-target % inhibition breakdown.

    Document counts are computed separately for literature and deposited
    records. Mixing them (dividing literature records by *all* documents)
    understates the per-paper figure, because deposited datasets carry many
    document ids of their own.
    """
    sql = f"""
        SELECT td.chembl_id                                           AS target,
               COUNT(*)                                               AS inh_total,
               SUM(CASE WHEN {src} =  ? THEN 1 ELSE 0 END)            AS inh_lit,
               SUM(CASE WHEN {src} <> ? THEN 1 ELSE 0 END)            AS inh_dep,
               COUNT(DISTINCT CASE WHEN {src} =  ? THEN act.doc_id END)   AS lit_docs,
               COUNT(DISTINCT CASE WHEN {src} <> ? THEN act.doc_id END)   AS dep_docs,
               COUNT(DISTINCT CASE WHEN {src} =  ? THEN act.molregno END) AS lit_mols,
               COUNT(DISTINCT act.molregno)                           AS all_mols,
               SUM(CASE WHEN d.year >= ? THEN 1 ELSE 0 END)           AS inh_recent,
               SUM(CASE WHEN d.year >= ? AND {src} = ? THEN 1 ELSE 0 END) AS inh_recent_lit,
               SUM(CASE WHEN act.standard_value IS NULL THEN 1 ELSE 0 END)  AS missing_value,
               SUM(CASE WHEN act.standard_value < 0 THEN 1 ELSE 0 END)      AS negative_value,
               SUM(CASE WHEN act.standard_value < 0 AND {src} =  ? THEN 1 ELSE 0 END) AS negative_lit,
               SUM(CASE WHEN act.standard_value < 0 AND {src} <> ? THEN 1 ELSE 0 END) AS negative_dep,
               SUM(CASE WHEN act.standard_relation IS NOT NULL
                         AND act.standard_relation <> '=' THEN 1 ELSE 0 END) AS censored_value
        {BASE_JOIN}
        WHERE td.chembl_id IN ({','.join('?' * len(ids))})
          AND act.standard_type = 'Inhibition'
        GROUP BY td.chembl_id
    """
    L = SRC_LITERATURE
    return con.execute(sql, [L, L, L, L, L, year, year, L, L, L, *ids]).fetchall()


def totals(con, ids, year):
    sql = f"""
        SELECT td.chembl_id                                    AS target,
               COUNT(*)                                        AS total,
               SUM(CASE WHEN d.year >= ? THEN 1 ELSE 0 END)    AS total_recent
        {BASE_JOIN}
        WHERE td.chembl_id IN ({','.join('?' * len(ids))})
        GROUP BY td.chembl_id
    """
    return con.execute(sql, [year, *ids]).fetchall()


def inhibition_documents(con, ids, year, src):
    """One row per (target, document) for literature-sourced % inhibition.

    This is the candidate list for the supplementary-information sampling:
    each row is a real paper, with its DOI, that contributed % inhibition for
    one of the targets.
    """
    sql = f"""
        SELECT td.chembl_id                     AS target,
               d.chembl_id                      AS doc,
               d.doi, d.pubmed_id, d.year, d.journal,
               SUBSTR(COALESCE(d.title, ''), 1, 200) AS title,
               COUNT(*)                         AS n_records,
               COUNT(DISTINCT act.molregno)     AS n_molecules
        {BASE_JOIN}
        WHERE td.chembl_id IN ({','.join('?' * len(ids))})
          AND act.standard_type = 'Inhibition'
          AND {src} = ?
          AND d.year >= ?
        GROUP BY td.chembl_id, d.chembl_id
        ORDER BY n_records DESC
    """
    return con.execute(sql, [*ids, SRC_LITERATURE, year]).fetchall()


def inhibition_assays(con, ids, year, src):
    """Distinct assay descriptions behind literature % inhibition records.

    The test concentration and, where stated, the panel name live in this free
    text. Reading the distribution is the cheapest way to see what the metadata
    actually looks like before writing any parser for it.
    """
    sql = f"""
        SELECT td.chembl_id                               AS target,
               SUBSTR(COALESCE(a.description, ''), 1, 250) AS description,
               COUNT(*)                                   AS n_records,
               COUNT(DISTINCT act.doc_id)                 AS n_documents
        {BASE_JOIN}
        WHERE td.chembl_id IN ({','.join('?' * len(ids))})
          AND act.standard_type = 'Inhibition'
          AND {src} = ?
          AND d.year >= ?
        GROUP BY td.chembl_id, description
        ORDER BY n_records DESC
    """
    return con.execute(sql, [*ids, SRC_LITERATURE, year]).fetchall()


# --------------------------------------------------------------------------- #
def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {path}  ({len(rows)} rows)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="path to chembl_XX.db")
    ap.add_argument("--year", type=int, default=2019,
                    help="start of the publication-year window (default 2019)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--schema", action="store_true",
                    help="print table columns and exit")
    args = ap.parse_args()

    if not os.path.exists(args.db):
        sys.exit(f"Database not found: {args.db}")

    con = sqlite3.connect(args.db)
    con.execute("PRAGMA query_only = ON")

    if args.schema:
        show_schema(con)
        return 0

    print("Resolving schema...")
    src = locate_src(con)

    label = {cid: lab for cid, lab, _ in TARGETS}
    gene = {cid: g for cid, _, g in TARGETS}
    ids = [cid for cid, _, _ in TARGETS]

    print(f"\nQuerying {len(ids)} targets, year window from {args.year}...")

    tot = {r[0]: (r[1], r[2]) for r in totals(con, ids, args.year)}

    # ---- by standard_type ----
    rows = []
    for target, std_type, n_all, n_recent in by_type(con, ids, args.year):
        rows.append([label.get(target, target), gene.get(target, ""), target,
                     std_type, n_all, n_recent])
    rows.sort(key=lambda r: (r[0], -r[4]))
    write_csv(os.path.join(args.outdir, "survey_by_type.csv"),
              ["label", "gene", "chembl_id", "standard_type",
               "n_all_years", f"n_from_{args.year}"], rows)

    # ---- inhibition summary ----
    header = ["label", "gene", "chembl_id", "total_activities",
              f"total_from_{args.year}",
              "inhibition", f"inhibition_from_{args.year}",
              "inh_literature", "inh_deposited",
              f"inh_literature_from_{args.year}",
              "literature_documents", "deposited_documents",
              "literature_records_per_document",
              "literature_molecules", "all_molecules",
              "missing_value", "negative_value",
              "negative_literature", "negative_deposited", "censored_value",
              "inhibition_fraction", "literature_inhibition_fraction"]
    rows = []
    for r in inhibition_summary(con, ids, args.year, src):
        (target, inh, lit, dep, ldocs, ddocs, lmols, amols,
         recent, recent_lit, missing, negative, neg_lit, neg_dep, censored) = r
        t_all, t_recent = tot.get(target, (0, 0))
        rows.append([
            label.get(target, target), gene.get(target, ""), target,
            t_all, t_recent, inh, recent, lit, dep, recent_lit,
            ldocs, ddocs,
            round(lit / ldocs, 2) if ldocs else "",
            lmols, amols, missing, negative, neg_lit, neg_dep, censored,
            round(inh / t_all, 4) if t_all else "",
            round(lit / t_all, 4) if t_all else "",
        ])
    rows.sort(key=lambda r: -(r[5] or 0))
    write_csv(os.path.join(args.outdir, "survey_inhibition.csv"), header, rows)

    # ---- per-document candidates ----
    rows = [[label.get(r[0], r[0]), *r[1:]]
            for r in inhibition_documents(con, ids, args.year, src)]
    write_csv(os.path.join(args.outdir, "inhibition_documents.csv"),
              ["target", "doc_chembl_id", "doi", "pubmed_id", "year",
               "journal", "title", "n_records", "n_molecules"], rows)

    # ---- assay description distribution ----
    rows = [[label.get(r[0], r[0]), *r[1:]]
            for r in inhibition_assays(con, ids, args.year, src)]
    write_csv(os.path.join(args.outdir, "inhibition_assays.csv"),
              ["target", "assay_description", "n_records", "n_documents"], rows)

    con.close()
    print("\nDone.")
    print("Next: inhibition_documents.csv lists real papers with DOIs that "
          "contributed % inhibition. Sample from it to measure how much of "
          "each paper's panel ChEMBL actually captured.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
