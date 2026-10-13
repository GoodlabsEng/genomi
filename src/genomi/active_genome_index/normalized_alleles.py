"""Transactional derived alleles inside an AGI, linked to unchanged observations.

Normalization is supplied by the preparation pipeline with its reference identity.
This module owns the index layout and readiness contract; the original records,
including their phase, confidence, ALT order and coordinates, remain authoritative.
"""
from __future__ import annotations

import json
import re
import sqlite3
from collections.abc import Iterable, Mapping
from typing import Any
from .array_genotypes import array_snp_call

NORMALIZED_ALLELE_SCHEMA_VERSION = 1
_SHA256 = re.compile(r"[a-f0-9]{64}")
_COLUMNS = ("batch_id", "chrom", "pos", "ref", "alt", "record_rowid",
            "original_alt", "original_allele_index", "match_mode", "kind")


def normalized_allele_index_metadata(connection: sqlite3.Connection) -> dict[str, Any]:
    """Return completed enrichment metadata, refusing an incomplete or old layout."""
    row = connection.execute("select value from metadata where key='normalized_allele_index'").fetchone()
    value = json.loads(row[0]) if row else None
    fields = {"schemaVersion", "complete", "genomeBuild", "normalizationVersion",
              "referenceSha256", "sourceSha256", "observedAlleleCount", "excludedContigs"}
    if (not isinstance(value, dict) or set(value) != fields
            or type(value.get("schemaVersion")) is not int or value["schemaVersion"] != NORMALIZED_ALLELE_SCHEMA_VERSION
            or value.get("complete") is not True or value.get("genomeBuild") not in {"GRCh37", "GRCh38"}
            or not isinstance(value.get("normalizationVersion"), str) or not 1 <= len(value["normalizationVersion"]) <= 100
            or not isinstance(value.get("referenceSha256"), str) or not _SHA256.fullmatch(value["referenceSha256"])
            or not isinstance(value.get("sourceSha256"), str) or not _SHA256.fullmatch(value["sourceSha256"])
            or type(value.get("observedAlleleCount")) is not int or not 0 <= value["observedAlleleCount"] <= 2**53-1
            or not isinstance(value.get("excludedContigs"), list) or len(value["excludedContigs"]) > 3):
        raise ValueError("normalized_allele_index_incomplete")
    seen = set()
    for exclusion in value["excludedContigs"]:
        if (not isinstance(exclusion, dict) or set(exclusion) != {"reason", "contigCount", "indexedRecordCount", "observedAlleleCount"}
                or exclusion.get("reason") not in {"supplemental_decoy", "supplemental_hla", "supplemental_viral"}
                or exclusion["reason"] in seen
                or any(type(exclusion[k]) is not int or not 0 <= exclusion[k] <= 2**53-1
                       for k in ("contigCount", "indexedRecordCount", "observedAlleleCount"))
                or exclusion["contigCount"] < 1 or exclusion["indexedRecordCount"] < exclusion["contigCount"]):
            raise ValueError("normalized_allele_index_incomplete")
        seen.add(exclusion["reason"])
    for name in ("normalized_alleles", "chromosome_aliases"):
        kind = connection.execute("select type,rootpage from sqlite_schema where name=?", (name,)).fetchone()
        if kind is None or kind[0] != "table" or kind[1] <= 0:
            raise ValueError("normalized_allele_index_incomplete")
    names = {row[0] for row in connection.execute("select name from sqlite_schema where type='index'")}
    if not {"normalized_seq_idx", "normalized_array_idx", "chromosome_aliases_canonical_idx"} <= names:
        raise ValueError("normalized_allele_index_incomplete")
    return value


def write_normalized_allele_index(
    connection: sqlite3.Connection,
    *,
    expected_aliases: Mapping[str, str],
    normalized_rows: Iterable[tuple],
    metadata: Mapping[str, Any],
) -> dict[str, Any]:
    """Install one completed derived index atomically; never rewrite source rows.

    The caller may normalize into temporary tables on this connection before calling.
    A savepoint keeps interrupted or invalid enrichment from marking the AGI ready.
    Existing completed enrichments cannot be overwritten with different source bindings.
    """
    if connection.execute("select 1 from metadata where key='normalized_allele_index'").fetchone():
        raise ValueError("normalized_allele_index_already_exists")
    if (metadata.get("schemaVersion") != NORMALIZED_ALLELE_SCHEMA_VERSION
            or metadata.get("genomeBuild") not in {"GRCh37", "GRCh38"}
            or not isinstance(metadata.get("referenceSha256"), str) or not _SHA256.fullmatch(metadata["referenceSha256"])
            or not isinstance(metadata.get("sourceSha256"), str) or not _SHA256.fullmatch(metadata["sourceSha256"])
            or not isinstance(metadata.get("normalizationVersion"), str)
            or not isinstance(metadata.get("excludedContigs"), list)):
        raise ValueError("invalid_normalized_allele_metadata")
    expected_aliases = dict(expected_aliases)
    if not expected_aliases or any(not k or not v for k, v in expected_aliases.items()):
        raise ValueError("invalid_normalized_allele_aliases")
    connection.create_function("array_snp_call", 2, array_snp_call, deterministic=True)
    connection.execute("savepoint normalized_allele_write")
    try:
        connection.execute("create table chromosome_aliases(alias text primary key, canonical text not null)")
        connection.execute("create index chromosome_aliases_canonical_idx on chromosome_aliases(canonical)")
        connection.execute("""create table normalized_alleles(
            batch_id text primary key, chrom text not null, pos integer not null,
            ref text not null, alt text not null, record_rowid integer not null,
            original_alt text not null, original_allele_index integer,
            match_mode text not null, kind text not null
        )""")
        connection.execute("create index normalized_seq_idx on normalized_alleles(chrom,pos,ref,alt) where kind='variant_call'")
        connection.execute("create index normalized_array_idx on normalized_alleles(chrom,pos,alt) where kind='array_call'")
        connection.executemany("insert into chromosome_aliases values (?,?)", expected_aliases.items())
        count = 0
        canonical_contigs = set(expected_aliases.values())
        pending = []
        for row in normalized_rows:
            if (len(row) != len(_COLUMNS) or type(row[2]) is not int or row[2] < 1
                    or row[1] not in canonical_contigs or type(row[5]) is not int or row[5] < 1
                    or row[9] not in {"variant_call", "array_call"}
                    or (row[9] == "variant_call" and (type(row[7]) is not int or row[7] < 1))
                    or (row[9] == "array_call" and row[7] is not None)):
                raise ValueError("invalid_normalized_allele_row")
            pending.append(row)
            count += 1
            if len(pending) == 10000:
                connection.executemany("insert into normalized_alleles values(?,?,?,?,?,?,?,?,?,?)", pending)
                pending.clear()
        if pending:
            connection.executemany("insert into normalized_alleles values(?,?,?,?,?,?,?,?,?,?)", pending)
        if connection.execute("""select 1 from normalized_alleles a left join records r
                on r.rowid=a.record_rowid where r.rowid is null or a.kind != r.record_kind limit 1""").fetchone():
            raise ValueError("normalized_allele_source_missing")
        if "observedAlleleCount" in metadata and count != metadata["observedAlleleCount"]:
            raise ValueError("normalized_allele_count_mismatch")
        invalid = connection.execute("""select 1 from normalized_alleles a join records r on r.rowid=a.record_rowid
            left join chromosome_aliases ca on ca.alias=r.chrom
            where ca.canonical is null or ca.canonical!=a.chrom
            or not exists (select 1 from json_each(r.observed_alleles) o where upper(o.value)=a.original_alt)
            or (a.kind='variant_call' and not exists(
                select 1 from json_each('["'||replace(upper(r.alt),',','","')||'"]') alt
                where alt.key=a.original_allele_index-1 and alt.value=a.original_alt
                      and a.original_alt!=upper(r.ref)))
            or (a.kind='array_call' and (a.original_allele_index is not null
                or not array_snp_call(r.genotype,r.observed_alleles)
                or a.original_alt not in ('A','C','G','T') or a.pos!=r.pos or a.alt!=a.original_alt))
            limit 1""").fetchone()
        if invalid:
            raise ValueError("normalized_allele_source_allele_mismatch")
        value = {**metadata, "observedAlleleCount": count, "complete": True}
        connection.execute("insert into metadata(key,value) values('normalized_allele_index',?)",
                           (json.dumps(value, sort_keys=True),))
        normalized_allele_index_metadata(connection)
        connection.execute("release normalized_allele_write")
        return value
    except BaseException:
        connection.execute("rollback to normalized_allele_write")
        connection.execute("release normalized_allele_write")
        raise
