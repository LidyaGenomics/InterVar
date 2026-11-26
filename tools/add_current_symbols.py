#!/usr/bin/env python3
"""Add current HGNC symbols to InterVar's gene-keyed lists (in place).

InterVar's intervardb gene lists were frozen in 2021 and are matched by gene
symbol against the gene database's symbols. Genes HGNC renamed since then
(ARNTL -> BMAL1, DDX58 -> RIGI, AARS -> AARS1, ...) silently stop matching, so
e.g. a loss-of-function variant in such a gene never gets PVS1. For every list
entry that is an HGNC previous symbol of exactly one approved gene, this adds a
line with the approved symbol and keeps the old line, so both spellings match.

Lists: PVS1.LOF.genes, PP2.genes, BP1.genes (one symbol per line) and
PM1_domains_with_benigns (gene in column 2), for hg19 and hg38.

    python3 tools/add_current_symbols.py hgnc_complete_set.txt intervardb/
"""
import csv
import os
import sys
from collections import defaultdict

ONE_PER_LINE = ["PVS1.LOF.genes", "PP2.genes", "BP1.genes"]
TABLES = {"PM1_domains_with_benigns": 1}  # file -> 0-based gene column


def previous_to_current(hgnc_path):
    approved = set()
    prev = defaultdict(set)
    alias = defaultdict(set)
    with open(hgnc_path, encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            if rec.get("status") != "Approved":
                continue
            approved.add(rec["symbol"])
            for col, into in (("prev_symbol", prev), ("alias_symbol", alias)):
                for p in (rec.get(col) or "").strip('"').split("|"):
                    if p:
                        into[p].add(rec["symbol"])
    # Only unambiguous renames of symbols that are no longer approved themselves;
    # a previous symbol wins, an alias is used only when it is not a previous
    # symbol of anything (HDGFRP3 -> HDGFL3, ZAK -> MAP3K20).
    out = {p: next(iter(c)) for p, c in alias.items() if len(c) == 1 and p not in approved}
    out.update({p: next(iter(c)) for p, c in prev.items() if len(c) == 1 and p not in approved})
    return out


def update_list(path, mapping):
    with open(path, encoding="utf-8", errors="surrogateescape") as fh:
        lines = [l.rstrip("\n") for l in fh]
    present = {l.strip() for l in lines}
    added = []
    for l in lines:
        cur = mapping.get(l.strip())
        if cur and cur not in present:
            present.add(cur)
            added.append(cur)
    if added:
        with open(path, "a", encoding="utf-8") as fh:
            fh.write("".join(f"{s}\n" for s in added))
    return added


def update_table(path, col, mapping):
    with open(path, encoding="utf-8", errors="surrogateescape") as fh:
        lines = [l.rstrip("\n") for l in fh]
    present = {tuple(l.split("\t")) for l in lines}
    added = []
    for l in lines[1:]:
        t = l.split("\t")
        if len(t) > col and t[col] in mapping:
            n = t[:col] + [mapping[t[col]]] + t[col + 1 :]
            if tuple(n) not in present:
                present.add(tuple(n))
                added.append("\t".join(n))
    if added:
        with open(path, "a", encoding="utf-8", errors="surrogateescape") as fh:
            fh.write("".join(f"{l}\n" for l in added))
    return added


def main(hgnc_path, db_dir):
    mapping = previous_to_current(hgnc_path)
    for build in ("hg19", "hg38"):
        for name in ONE_PER_LINE:
            p = os.path.join(db_dir, f"{name}.{build}")
            if os.path.exists(p):
                print(f"{name}.{build}: +{len(update_list(p, mapping))} current symbols")
        for name, col in TABLES.items():
            p = os.path.join(db_dir, f"{name}.{build}")
            if os.path.exists(p):
                print(f"{name}.{build}: +{len(update_table(p, col, mapping))} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
