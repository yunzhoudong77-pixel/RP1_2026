#!/usr/bin/env python3
"""Search Europe PMC for publications likely to contain assay-panel data.

Why not just the vendor names
-----------------------------
Two panel-bearing papers inspected by hand (J Med Chem 2019, 62, 2750 and
J Med Chem 2020, jm0c00814) name no vendor anywhere in the text: no
"Eurofins", no "Cerep", no "SafetyScreen44". Brand-name queries therefore
undercount badly. The query set below is organised in tiers so that the
precision of each tier can be measured separately once papers are checked:

    vendor      product and company names          high precision, low recall
    phrasing    the wording panel reports use      high precision
    descriptive generic panel language             medium
    readout     readout + concentration combos     low precision, high recall
    transporter regulatory transporter panels      separate data stream

Two filters are deliberately NOT applied:

* HAS_SUPPL is not required. J Med Chem 2019, 62, 2750 carries a full
  off-target panel as a main-text table and has no supplementary files at all;
  requiring supplementary material would discard exactly the papers we want.
* OPEN_ACCESS is not required either, but is recorded, since it determines
  whether supplementary files can be fetched programmatically.

Usage
-----
    pip install requests
    python epmc_search.py                      # all tiers, 2019-2026
    python epmc_search.py --tier vendor phrasing
    python epmc_search.py --counts-only        # just hit counts per query
    python epmc_search.py --from 2015 --to 2026

Outputs
-------
    epmc_query_counts.csv    one row per query: hit count, how many are OA,
                             how many carry supplementary files
    epmc_candidates.csv      one row per publication, deduplicated, with the
                             list of queries that matched it
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

try:
    import requests
except ImportError:
    sys.exit("Missing dependency. Run: pip install requests")

REST = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CACHE = os.path.join(os.environ.get("LOCALAPPDATA", "."), "epmc_search", "epmc_cache.json")
os.makedirs(os.path.dirname(CACHE), exist_ok=True)


def cache_save() -> None:
    tmp = CACHE + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(_cache, fh)
    for attempt in range(10):          # 万一仍被占用，等一下再试
        try:
            os.replace(tmp, CACHE)
            return
        except PermissionError:
            time.sleep(0.5 * (attempt + 1))
    print("    warning: cache not saved this round (file locked)")
PAGE = 1000
PAUSE = 0.34
RETRIES = 4

# --------------------------------------------------------------------------- #
# Query set
# --------------------------------------------------------------------------- #
QUERIES: dict[str, list[str]] = {

    # Tightest set: vendor names that cannot mean anything else, plus the exact
    # wording Cerep-style panel reports use. Run this tier alone for a pilot --
    # it returns ~800 publications for 2019-2026, ~535 of them open access with
    # supplementary files.
    "core": [
        '"inhibition of control specific binding"',
        '"control specific binding"',
        '"percent inhibition of control"',
        '"SafetyScreen44"',
        '"Safety Screen 44"',
        '"SafetyScreen 44"',
        '"Eurofins Cerep"',
        '"Eurofins Panlabs"',
        '"off-target panel"',
        '"panel of off-targets"',
        '"safety pharmacology panel"',
        '"receptor binding panel"',
        '"GPCRScan"',
    ],

    # Vendor and product names that are specific but also appear in passing,
    # e.g. as the supplier of a single assay rather than a panel.
    "vendor": [
        '"Cerep"',
        '"PDSP"',
        '"Psychoactive Drug Screening Program"',
    ],

    # Generic panel language, for papers that name no vendor at all -- which
    # both hand-checked panel papers turned out to do.
    "descriptive": [
        '"off-target profiling"',
        '"selectivity panel"',
        '"secondary pharmacology"',
        '"binding panel"',
        '"counterscreen"',
        '"broad selectivity"',
    ],

    # The FDA / EMA / ICH M12 mandatory transporter panel: a separate stream.
    "transporter": [
        '"OATP1B1" AND "OAT1" AND "OCT2" AND "MATE1"',
        '"OATP1B1" AND "OATP1B3" AND "BCRP" AND "P-gp"',
        '"MATE2-K" AND "inhibition"',
        '"ICH M12"',
        '"transporter inhibition" AND "panel"',
        '"Solvo Biotechnology"',
        '"Optivia"',
    ],

    # Deliberately kept out of the default run; enable with
    # --tier wide to measure how much the tighter tiers miss.
    "wide": [
        '"percent of control"',
        '"single concentration" AND ("off-target" OR "selectivity")',
        '"percent inhibition" AND "selectivity" AND "panel"',
        '"counter-screen"',
    ],
}

# Queries removed after measuring them, and why. Kept as a record so the same
# mistakes are not repeated, and because the first of them is a property of the
# search engine rather than of the literature.
#
#   "% of control"                 5,088,974 hits. Europe PMC does not treat
#   "% inhibition of control"      "%" as part of a phrase, so the query
#                                  degrades to a search for "of control".
#                                  Always write "percent", never "%".
#
#   "BioPrint"                     1,455 hits, of which 1,441 also mention 3D
#                                  printing, bioprinting or tissue: these are
#                                  3D bioprinting papers, not Cerep BioPrint.
#
#   "broad panel"                  4,598 hits, 2,895 of them matched by this
#                                  query alone. "A broad panel of experts",
#                                  "a broad panel of cell lines", and so on.
#
#   "residual activity" AND "panel"        3,746 hits, 2,848 exclusive.
#   "single-point" AND "screening" ...     4,767 hits, 2,653 exclusive.
#   "single concentration" AND "inhibition" 3,604 hits, 2,514 exclusive.
#                                  All three are too loose to be useful without
#                                  a panel term alongside them; narrowed
#                                  versions live in the "wide" tier.


# --------------------------------------------------------------------------- #
# HTTP with cache and retry
# --------------------------------------------------------------------------- #
_cache: dict = {}
_session = requests.Session()
_session.headers["User-Agent"] = "avoidome-sdp-survey/1.0"


def cache_load() -> None:
    global _cache
    if os.path.exists(CACHE):
        try:
            _cache = json.load(open(CACHE, encoding="utf-8"))
            print(f"loaded {len(_cache)} cached responses")
        except Exception:
            _cache = {}


def cache_save() -> None:
    tmp = CACHE + ".tmp"
    json.dump(_cache, open(tmp, "w", encoding="utf-8"))
    os.replace(tmp, CACHE)


def fetch(params: dict) -> dict | None:
    key = json.dumps(params, sort_keys=True)
    if key in _cache:
        return _cache[key]
    for attempt in range(1, RETRIES + 1):
        try:
            r = _session.get(REST, params=params, timeout=90)
            if r.status_code == 429:
                time.sleep(5 * attempt)
                continue
            r.raise_for_status()
            payload = r.json()
            _cache[key] = payload
            time.sleep(PAUSE)
            return payload
        except Exception as exc:
            print(f"    retry {attempt}/{RETRIES}: {type(exc).__name__}")
            time.sleep(2 * attempt)
    print(f"    FAILED: {params.get('query','')[:70]}")
    return None


# --------------------------------------------------------------------------- #
def window(q: str, y_from: int, y_to: int) -> str:
    return f"({q}) AND (FIRST_PDATE:[{y_from}-01-01 TO {y_to}-12-31])"


def count(q: str) -> int | None:
    p = fetch({"query": q, "format": "json", "pageSize": 1})
    return p["hitCount"] if p else None


def harvest(q: str, cap: int) -> list[dict]:
    """Page through a query with cursorMark, up to `cap` records."""
    out, cursor = [], "*"
    while len(out) < cap:
        p = fetch({"query": q, "format": "json", "resultType": "core",
                   "pageSize": min(PAGE, cap - len(out)), "cursorMark": cursor})
        if not p:
            break
        batch = p.get("resultList", {}).get("result", [])
        if not batch:
            break
        out.extend(batch)
        nxt = p.get("nextCursorMark")
        if not nxt or nxt == cursor:
            break
        cursor = nxt
    return out


def row_of(rec: dict) -> dict:
    j = rec.get("journalInfo") or {}
    return {
        "pmid": rec.get("pmid", ""),
        "pmcid": rec.get("pmcid", ""),
        "doi": rec.get("doi", ""),
        "year": rec.get("pubYear", ""),
        "journal": ((j.get("journal") or {}).get("title") or "")[:60],
        "title": " ".join((rec.get("title") or "").split())[:240],
        "open_access": rec.get("isOpenAccess", ""),
        "has_suppl": rec.get("hasSuppl", ""),
        "source": rec.get("source", ""),
        "cited_by": rec.get("citedByCount", ""),
    }


def key_of(rec: dict) -> str:
    return (rec.get("doi") or rec.get("pmid")
            or rec.get("pmcid") or rec.get("id") or repr(rec))[:120]


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", nargs="*",
                    default=["core", "vendor", "descriptive", "transporter"],
                    choices=list(QUERIES), help="which query tiers to run")
    ap.add_argument("--from", dest="y_from", type=int, default=2019)
    ap.add_argument("--to", dest="y_to", type=int, default=2026)
    ap.add_argument("--cap", type=int, default=3000,
                    help="max records fetched per query (default 3000)")
    ap.add_argument("--counts-only", action="store_true",
                    help="report hit counts without downloading records")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    cache_load()

    jobs = [(tier, q) for tier in args.tier for q in QUERIES[tier]]
    print(f"{len(jobs)} queries, window {args.y_from}-{args.y_to}\n")

    # ---- counts -------------------------------------------------------- #
    counts = []
    for tier, q in jobs:
        wq = window(q, args.y_from, args.y_to)
        n = count(wq)
        n_oa = count(wq + " AND OPEN_ACCESS:y")
        n_sup = count(wq + " AND HAS_SUPPL:y")
        counts.append({"tier": tier, "query": q, "hits": n,
                       "open_access": n_oa, "has_suppl": n_sup})
        print(f"  [{tier:<11}] {str(n):>6} hits "
              f"({n_oa} OA, {n_sup} with SI)   {q[:60]}")
        cache_save()

    path = os.path.join(args.outdir, "epmc_query_counts.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["tier", "query", "hits",
                                           "open_access", "has_suppl"])
        w.writeheader()
        w.writerows(counts)
    print(f"\nwrote {path}")

    if args.counts_only:
        return 0

    # ---- records ------------------------------------------------------- #
    print("\nfetching records...")
    merged: dict[str, dict] = {}
    for tier, q in jobs:
        wq = window(q, args.y_from, args.y_to)
        recs = harvest(wq, args.cap)
        print(f"  [{tier:<11}] {len(recs):>5} records   {q[:60]}")
        for rec in recs:
            k = key_of(rec)
            if k not in merged:
                merged[k] = row_of(rec)
                merged[k]["matched"] = []
                merged[k]["tiers"] = set()
            merged[k]["matched"].append(q)
            merged[k]["tiers"].add(tier)
        cache_save()

    rows = []
    for v in merged.values():
        v["n_queries_matched"] = len(v["matched"])
        v["matched"] = " | ".join(v["matched"])[:400]
        v["tiers"] = ",".join(sorted(v["tiers"]))
        rows.append(v)
    # papers matched by several queries first: most likely to be real panels
    rows.sort(key=lambda r: (-r["n_queries_matched"], r["year"]), reverse=False)
    rows.sort(key=lambda r: -r["n_queries_matched"])

    cols = ["doi", "pmid", "pmcid", "year", "journal", "title",
            "open_access", "has_suppl", "source", "cited_by",
            "tiers", "n_queries_matched", "matched"]
    path = os.path.join(args.outdir, "epmc_candidates.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    oa_si = sum(1 for r in rows if r["open_access"] == "Y" and r["has_suppl"] == "Y")
    print(f"\nwrote {path}")
    print(f"  {len(rows)} unique publications")
    print(f"  {oa_si} open access with supplementary files (fetchable directly)")
    print(f"  {sum(1 for r in rows if r['n_queries_matched'] >= 3)} matched by "
          f"3 or more queries (strongest candidates)")
    cache_save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
