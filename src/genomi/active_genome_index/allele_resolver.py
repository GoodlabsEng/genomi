"""Batch resolution against normalized AGI alleles and original observations.

Resolve public targets against one enriched, immutable Active Genome Index. Normalized
entries establish complete sequence-allele identity; original calls establish genotype,
dosage, phase and confidence. A projected base never establishes an exact allele match.

The caller proves reference identity against normalized_allele_index.referenceSha256
before entering this API. This module checks reference anchors without hashing an entire
FASTA on each query. The reader wrapper must apply its usual authorization/readiness gate.

Normalized index metadata used here: schemaVersion=1, complete=True, genomeBuild,
normalizationVersion, referenceSha256 and observedAlleleCount. The existing parse schema
is unchanged. The internal exact-hit ordinals describe the complete candidate inventory;
only the first 100 relevant targets receive source assessment. The caller orders public
targets by complete clinical priority before resolving and presents explicit omissions.


"""
from __future__ import annotations

from collections import defaultdict
from contextlib import ExitStack
import json
import math
from pathlib import Path
import re
import sqlite3
import time
from typing import TYPE_CHECKING

from ._agi_query import _index_record_to_dict
from .observations import observed_alleles_from_vcf_genotype
from .array_genotypes import classify_array_genotype
from .vcf import parse_record_fields

if TYPE_CHECKING:
    from .reader import ActiveGenomeIndexReader

MAX_SOURCE_RECORDS = 64
MAX_RESOLVED_TARGETS = 100
MAX_INPUT_TARGETS = 100_000
MAX_SOURCE_LOOKUP_BYTES = 8 * 1024 * 1024
MIN_DEPTH = 10
MIN_GENOTYPE_QUALITY = 20
GENOTYPE_PATTERN = re.compile(r"(?:[0-9]+|\.)(?:[/|](?:[0-9]+|\.))*\Z")
SEQUENCE_PATTERN = re.compile(r"[ACGTN]+\Z")
CALLED_SEQUENCE_PATTERN = re.compile(r"[ACGT]+\Z")
ARRAY_FORMATS = frozenset(("23andme", "nucleus"))


class AlleleResolutionUnavailable(ValueError):
    """Only an allowlisted reason is exposed; source or SQLite messages stay private."""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


def remaining_time(deadline: float) -> float:
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("allele_resolution_deadline_exceeded")
    return remaining


def _metadata(connection: sqlite3.Connection) -> dict:
    values = {}
    for key, value in connection.execute("select key,value from metadata"):
        try:
            values[key] = json.loads(value)
        except (TypeError, json.JSONDecodeError):
            values[key] = value
    return values


def _index_identity(metadata: dict, build: str | None, expected_reference_sha256: str) -> dict:
    identity = metadata.get("normalized_allele_index")
    if (not isinstance(identity, dict) or identity.get("genomeBuild") != build
            or identity.get("referenceSha256") != expected_reference_sha256):
        raise AlleleResolutionUnavailable("normalized_allele_index_incompatible")
    return identity


def sample_call_quality(record: dict) -> dict:
    """Only explicit sample FORMAT metrics count as sequencing confidence."""
    sample = record.get("sample")
    sample = sample if isinstance(sample, dict) else {}
    metrics = {
        key: int(value) if isinstance(value, str) and re.fullmatch(r"[0-9]+", value) else None
        for key in ("DP", "GQ", "MIN_DP")
        for value in (sample.get(key),)
    }
    depth, quality, minimum = metrics["DP"], metrics["GQ"], metrics["MIN_DP"]
    block = record.get("end", record["pos"]) > record["pos"] and not record.get("is_variant")
    missing, failed = [], []
    if depth is None:
        missing.append("depth_missing")
    if quality is None:
        missing.append("genotype_quality_missing")
    if block and minimum is None:
        missing.append("minimum_depth_missing")
    effective_depth = depth
    if block:
        effective_depth = min(depth, minimum) if depth is not None and minimum is not None else None
    if effective_depth is not None and effective_depth < MIN_DEPTH:
        failed.append("depth_below_threshold")
    if quality is not None and quality < MIN_GENOTYPE_QUALITY:
        failed.append("genotype_quality_below_threshold")
    record_filter, sample_filter = record.get("filter"), sample.get("FT")
    if record_filter in (None, "", "."):
        missing.append("record_filter_not_applied")
    elif record_filter != "PASS":
        failed.append("record_filter_failed")
    if sample_filter not in (None, "", ".", "PASS"):
        failed.append("sample_filter_failed")
    return {
        "quality": "low" if failed else "unknown" if missing else "pass",
        "quality_reasons": failed + missing,
        "depth": depth,
        "genotype_quality": quality,
        "minimum_depth": minimum,
        "record_filter": record_filter,
        "sample_filter": sample_filter,
    }


def assess_call(record: dict, *, array: bool = False) -> dict:
    genotype = record.get("genotype")
    sample = record.get("sample")
    sample = sample if isinstance(sample, dict) else {}
    separators = set(re.findall(r"[/|]", genotype)) if isinstance(genotype, str) else set()
    separator = next(iter(separators)) if len(separators) == 1 else "mixed" if separators else None
    result = {
        "status": "no_call", "raw_genotype": genotype, "ploidy": None,
        "alleles": [], "quality": "unknown", "zygosity": "unknown",
        "phased": separator == "|", "phase_separator": separator,
        "phase_set": sample.get("PS"), "phase_id": sample.get("PID"),
        "phase_genotype": sample.get("PGT"),
    }
    if array:
        result.update(
            probe_id=record.get("id") or record.get("rsid"),
            quality_reasons=["array_sequencing_metrics_unavailable"],
            depth=None, genotype_quality=None, minimum_depth=None,
            record_filter=None, sample_filter=None,
        )
        call = classify_array_genotype(genotype)
        if call.status == "unresolved" or (call.status == "called" and not call.snp_eligible):
            return {**result, "status": "unresolved", "reason": "array_genotype_unsupported_for_snp"}
        if call.status == "no_call" or record.get("record_kind") != "array_call":
            return result
        called = list(call.alleles)
        return {
            **result, "status": "called", "ploidy": len(called), "alleles": called,
            "zygosity": "haploid" if len(called) == 1 else
            "homozygous" if len(set(called)) == 1 else "heterozygous",
        }
    result.update(sample_call_quality(record))
    if not isinstance(genotype, str) or not GENOTYPE_PATTERN.fullmatch(genotype):
        result.update(
            quality="unknown",
            quality_reasons=[*result["quality_reasons"], "genotype_missing"],
        )
        return result
    tokens = re.split(r"[/|]", genotype)
    result["ploidy"] = len(tokens)
    if "." in tokens:
        result.update(
            quality="unknown",
            quality_reasons=[*result["quality_reasons"], "genotype_incomplete"],
        )
        return result
    alleles = [str(record.get("ref") or "").upper()]
    alleles.extend(str(record.get("alt") or "").upper().split(","))
    indexes = [int(token) for token in tokens]
    if any(index >= len(alleles) or not CALLED_SEQUENCE_PATTERN.fullmatch(alleles[index])
           for index in indexes):
        return {**result, "status": "unresolved"}
    called = [alleles[index] for index in indexes]
    if len(indexes) == 1:
        zygosity = "haploid_reference" if indexes[0] == 0 else "haploid_alternate"
    else:
        zygosity = "homozygous_reference" if all(index == 0 for index in indexes) else (
            "homozygous_alternate" if len(set(called)) == 1 else "heterozygous"
        )
    return {**result, "status": "called", "alleles": called, "zygosity": zygosity}


def assess_overlapping_calls(records: list[dict], *, array: bool = False,
                             probe_id: str | None = None) -> dict:
    calls = [assess_call(record, array=array) for record in records]
    if not calls:
        return {"status": "unresolved", "reason": "source_record_unavailable"}
    compared = list(zip(records, calls))
    missing_probes = False
    missing = []
    missing_ids = set()
    if array:
        complete = [(record, call) for record, call in compared if call["status"] == "called"]
        missing = [(record, call) for record, call in compared if call["status"] == "no_call"]
        called_ids = {record.get("id") or record.get("rsid") for record, _ in complete}
        missing_ids = {record.get("id") or record.get("rsid") for record, _ in missing}
        if (complete and missing and all(call["status"] in ("called", "no_call") for call in calls)
                and all(isinstance(identifier, str) and identifier for identifier in called_ids | missing_ids)
                and not called_ids & missing_ids):
            compared = complete
            missing_probes = True
    compared_records = [record for record, _ in compared]
    compared_calls = [call for _, call in compared]
    if any(call["status"] != "called" for call in compared_calls):
        status = "no_call" if any(call["status"] == "no_call" for call in compared_calls) else "unresolved"
    elif len({(record["pos"], str(record["ref"]).upper()) for record in compared_records}) != 1:
        status = "unresolved"
    elif len({tuple(sorted(call["alleles"])) for call in compared_calls}) == 1:
        status = "called"
    elif array or len({frozenset(str(record.get("alt") or "").upper().split(","))
                       for record in compared_records}) == 1:
        status = "conflicting"
    else:
        status = "unresolved"
    if status == "called" and missing_probes and probe_id in missing_ids:
        status = "no_call"
    order = sorted(zip(records, calls), key=lambda pair: (
        (pair[0].get("id") or pair[0].get("rsid")) != probe_id,
        pair[1]["status"] != "called",
    )) if array else list(zip(records, calls))
    provenance = {
        "source_calls": [call for _, call in order[:4]], "source_call_count": len(calls),
        "source_calls_omitted_count": max(0, len(calls) - 4),
    }
    if status != "called":
        return {"status": status, **provenance}
    result = dict(compared_calls[0])
    result["quality"] = "low" if any(call["quality"] == "low" for call in compared_calls) else (
        "unknown" if any(call["quality"] == "unknown" for call in compared_calls) else "pass"
    )
    result["quality_reasons"] = sorted({reason for call in compared_calls for reason in call["quality_reasons"]})
    for field in ("depth", "genotype_quality", "minimum_depth"):
        values = [call[field] for call in compared_calls]
        result[field] = min(values) if all(value is not None for value in values) else None
    # Genotype equivalence never grants phase equivalence or confidence from one duplicate.
    for field in ("record_filter", "sample_filter", "phased", "phase_separator", "phase_set", "phase_id", "phase_genotype"):
        values = [call[field] for call in compared_calls]
        if any(value != values[0] for value in values):
            result[field] = None
    if len({tuple(call["alleles"]) for call in compared_calls if call.get("phased") is True}) > 1:
        result.update(phased=None, phase_separator=None, phase_set=None,
                      phase_id=None, phase_genotype=None, phase_reason="duplicate_phase_disagreement")
    if len(calls) > 1:
        result.update(provenance)
    if missing_probes:
        result["missing_probe_count"] = len(missing)
    return result


class _OriginalSource:
    """One bounded BGZF reader; offsets bind to the AGI-owned canonical source."""

    def __init__(self, reader: ActiveGenomeIndexReader, metadata: dict, stack: ExitStack):
        from pysam.libcbgzf import BGZFile

        root = Path(reader.agi_path).resolve(strict=True).parent
        source_raw = metadata.get("vcf_path")
        if not isinstance(source_raw, str):
            raise AlleleResolutionUnavailable("canonical_source_unavailable")
        raw = Path(source_raw)
        source = raw.resolve(strict=True)
        if raw.is_symlink() or not source.is_relative_to(root) or not source.is_file():
            raise AlleleResolutionUnavailable("canonical_source_binding_invalid")
        raw_index = Path(str(raw) + ".gzi")
        index = raw_index.resolve(strict=True)
        if raw_index.is_symlink() or not index.is_relative_to(root) or not index.is_file():
            raise AlleleResolutionUnavailable("canonical_source_binding_invalid")
        stat = source.stat()
        if stat.st_size != metadata.get("vcf_size_bytes") or stat.st_mtime_ns != metadata.get("vcf_mtime_ns"):
            raise AlleleResolutionUnavailable("canonical_source_changed")
        self.handle = stack.enter_context(BGZFile(str(source), "rb"))
        self.cache = {}

    def record(self, indexed: dict, *, aliases: set[str], position: int,
               deadline: float, budget: list[int]) -> dict:
        remaining_time(deadline)
        if indexed.get("sample_index") != 0:
            raise AlleleResolutionUnavailable("multiple_source_samples")
        start = indexed["offset"]
        self.handle.seek(start)
        first = True
        while budget[0] > 0:
            remaining_time(deadline)
            offset = self.handle.tell()
            raw = self.handle.readline()
            budget[0] -= len(raw)
            if not raw or budget[0] < 0:
                break
            fields = raw.decode("utf-8").rstrip("\r\n").split("\t")
            if len(fields) != 10:
                raise AlleleResolutionUnavailable("invalid_source_record")
            if offset not in self.cache:
                source = parse_record_fields(
                    fields, sample_name=indexed.get("sample_name"), sample_index=0,
                    offset=offset, line_length=len(raw),
                )
                record = source.to_dict()
                # Do not inherit parser inference from INFO/DP, AD or PL.
                sample = record["sample"]
                for field, key in (("depth", "DP"), ("genotype_quality", "GQ")):
                    value = sample.get(key)
                    record[field] = int(value) if isinstance(value, str) and re.fullmatch(r"[0-9]+", value) else None
                gt = str(record.get("genotype") or ".")
                incomplete = not GENOTYPE_PATTERN.fullmatch(gt) or "." in re.split(r"[/|]", gt)
                record["record_kind"] = "no_call" if incomplete else (
                    "variant_call" if record["is_variant"] else "reference_block"
                )
                record["observed_alleles"] = observed_alleles_from_vcf_genotype(
                    record["ref"], str(record.get("alt") or "."), record.get("genotype"),
                ) or []
                record["rsid"] = record.get("id")
            else:
                record = self.cache[offset]
            if first and (record["pos"] != indexed["pos"] or record["chrom"] != indexed["chrom"]
                          or record["ref"] != indexed["ref"]):
                raise AlleleResolutionUnavailable("source_offset_mismatch")
            first = False
            if record["chrom"] not in aliases or record["pos"] > position:
                break
            if record["end"] >= position:
                self.cache[offset] = record
                return {**record, "record_rowid": indexed["record_rowid"]}
        raise AlleleResolutionUnavailable("original_block_not_resolved")


def _covering_index_rows(connection: sqlite3.Connection, canonical: str, position: int) -> list[dict]:
    # Point records use records_region_idx; spanning records use spans_region_idx and
    # the offset/sample index. Never ask native chromosome guessing to widen this scope.
    rows = connection.execute("""
        select r.rowid as record_rowid, r.*
        from chromosome_aliases ca
        cross join records r on r.chrom=ca.alias
        where ca.canonical=? and r.pos=?
        union all
        select r.rowid as record_rowid, r.*
        from chromosome_aliases ca
        cross join spans s on s.chrom=ca.alias
        cross join records r on r.offset=s.offset and r.sample_index=s.sample_index
        where ca.canonical=? and s.pos<? and s.end>=?
        limit ?
    """, (canonical, position, canonical, position, position, MAX_SOURCE_RECORDS + 1)).fetchall()
    if len(rows) > MAX_SOURCE_RECORDS:
        raise AlleleResolutionUnavailable("covering_record_limit")
    return [dict(row) for row in rows]


def _records_at(connection: sqlite3.Connection, canonical: str, position: int, *,
                array: bool, source: _OriginalSource | None, aliases: set[str], deadline: float) -> list[dict]:
    remaining_time(deadline)
    rows = _covering_index_rows(connection, canonical, position)
    records = {}
    budget = [MAX_SOURCE_LOOKUP_BYTES]
    for row in rows:
        remaining_time(deadline)
        if array:
            # Array offsets refer to original text, so the parse-owned rows are authoritative.
            record = _index_record_to_dict(row)
            record["record_rowid"] = row["record_rowid"]
        else:
            if source is None:
                raise AlleleResolutionUnavailable("canonical_source_unavailable")
            record = source.record(row, aliases=aliases, position=position, deadline=deadline, budget=budget)
        records[(record["offset"], record.get("sample_index", 0))] = record
    return list(records.values())


def _reference_status(records: list[dict], assessment: dict, target: dict, *,
                      canonical: str, source_format: str, fasta) -> tuple[str, str | None]:
    if len(target["ref"]) != 1 or len(target["alt"]) != 1:
        return "unavailable", "unsupported_reference_projection"
    if len(records) != 1:
        return "unavailable", "overlapping_source_records"
    record = records[0]
    tokens = re.split(r"[/|]", str(record.get("genotype") or "."))
    if any(token != "0" for token in tokens):
        return "no_call", "source_reference_genotype_unresolved"
    if assessment.get("quality") != "pass":
        reasons = assessment.get("quality_reasons", [])
        return "unavailable", reasons[0] if reasons else "source_quality_insufficient"
    # All-reference GT supports the validated first REF base in either format.
    # Interior projection requires a genuine gVCF reference-confidence row.
    if target["pos"] > record["pos"]:
        if source_format != "gvcf":
            return "unavailable", "plain_vcf_reference_projection_unsupported"
        sample = record.get("sample", {})
        if (record.get("alts") not in (["<NON_REF>"], ["<*>"])
                or sample.get("LAA") not in (None, "", ".", "1")):
            return "unavailable", "unsupported_source_reference_span"
        length = sample.get("LEN")
        if length not in (None, "", "."):
            if not isinstance(length, str) or not re.fullmatch(r"[0-9]{1,10}", length) or int(length) < 1:
                return "unavailable", "invalid_source_reference_length"
            if target["pos"] >= record["pos"] + int(length):
                return "unavailable", "outside_source_reference_span"
    try:
        target_ref = fasta.fetch(canonical, target["pos"] - 1, target["pos"]).upper()
        source_ref = fasta.fetch(canonical, record["pos"] - 1,
                                 record["pos"] - 1 + len(record["ref"])).upper()
    except (OSError, ValueError, KeyError):
        return "unavailable", "reference_sequence_unavailable"
    if target_ref != target["ref"] or source_ref != record["ref"].upper():
        return "unavailable", "reference_sequence_mismatch"
    return "reference_confirmed", None


def _site_projection(records: list[dict], target: dict) -> dict | None:
    """Return informational equal-length base projection, never exact identity or dosage."""
    if len(records) != 1 or len(target["ref"]) != 1 or len(target["alt"]) != 1:
        return None
    record = records[0]
    ref = str(record.get("ref") or "").upper()
    alts = str(record.get("alt") or "").upper().split(",")
    offset = target["pos"] - record["pos"]
    if len(ref) <= 1 or not 0 <= offset < len(ref):
        return None
    gt = str(record.get("genotype") or ".")
    if not GENOTYPE_PATTERN.fullmatch(gt) or "." in re.split(r"[/|]", gt):
        return None
    alleles = [ref, *alts]
    indexes = [int(token) for token in re.split(r"[/|]", gt)]
    if any(index >= len(alleles) or len(alleles[index]) != len(ref)
           or not CALLED_SEQUENCE_PATTERN.fullmatch(alleles[index]) for index in indexes):
        return {"status": "unavailable", "reason": "unsupported_sequence_projection",
                "exact_sequence_allele": False}
    return {
        "status": "projected", "allele_bases": [alleles[index][offset] for index in indexes],
        "exact_sequence_allele": False,
        "limitation": "These bases come from a larger source change and do not identify the complete requested allele.",
    }


def _combined_assessment(records: list[dict], *, array: bool, probe_id: str | None,
                         normalized_source_loci: set[tuple]) -> dict:
    assessment = assess_overlapping_calls(records, array=array, probe_id=probe_id)
    if len(normalized_source_loci) <= 1:
        return assessment
    calls = [assess_call(record, array=array) for record in records]
    return {
        "status": "unresolved", "reason": "multiple_normalized_source_representations",
        "source_calls": calls[:4], "source_call_count": len(calls),
        "source_calls_omitted_count": max(0, len(calls) - 4),
        "source_loci": [dict(zip(("chrom", "pos", "ref"), locus))
                        for locus in sorted(normalized_source_loci)[:4]],
        "source_loci_omitted_count": max(0, len(normalized_source_loci) - 4),
    }


def _target_dosage(hits: list[dict], records: list[dict], assessment: dict, target: dict,
                   *, array: bool) -> int | None:
    if assessment.get("status") != "called":
        return None
    if array:
        alleles = assessment.get("alleles", [])
        return alleles.count(target["alt"]) if all(base in "ACGT" for base in alleles) else None
    indexes_by_row = defaultdict(set)
    for hit in hits:
        index = hit.get("original_allele_index")
        if type(index) is not int or index < 1:
            raise AlleleResolutionUnavailable("normalized_source_allele_invalid")
        indexes_by_row[hit["record_rowid"]].add(index)
    dosages = set()
    for record in records:
        if record["record_rowid"] not in indexes_by_row:
            continue
        indexes = indexes_by_row[record["record_rowid"]]
        alleles = str(record.get("alt") or "").upper().split(",")
        for hit in (item for item in hits if item["record_rowid"] == record["record_rowid"]):
            index = hit["original_allele_index"]
            if index > len(alleles) or alleles[index - 1] != str(hit["original_alt"]).upper():
                raise AlleleResolutionUnavailable("normalized_source_allele_invalid")
        tokens = re.split(r"[/|]", str(record.get("genotype") or "."))
        dosages.add(sum(int(token) in indexes for token in tokens))
    return next(iter(dosages)) if len(dosages) == 1 else None


def _resolve_target(connection: sqlite3.Connection, target: dict, hits: list[dict], *,
                    canonical: str | None, source_format: str, source: _OriginalSource | None,
                    fasta, deadline: float, cache: dict) -> dict:
    result = {
        "target": target, "exact_match": bool(hits), "target_dosage": None,
        "source_records": [], "call_assessment": {"status": "unresolved"},
        "status": "unavailable", "reference_reason": None,
    }
    if len(hits) > MAX_SOURCE_RECORDS:
        result["reference_reason"] = "source_call_limit"
        return result
    if canonical is None:
        result["reference_reason"] = "unsupported_target_contig"
        return result
    identity_anchor = target["ref"] == target["alt"] and bool(target.get("rsid"))
    array = source_format in ARRAY_FORMATS
    if array and not identity_anchor and (
            len(target["ref"]) != 1 or len(target["alt"]) != 1
            or target["ref"] not in "ACGT" or target["alt"] not in "ACGT"):
        result["exact_match"] = False
        result["reference_reason"] = "array_allele_model_unsupported"
        return result
    if identity_anchor:
        result["exact_match"] = False
    if not identity_anchor and (not SEQUENCE_PATTERN.fullmatch(target["ref"])
                                or not SEQUENCE_PATTERN.fullmatch(target["alt"])):
        result["reference_reason"] = "unsupported_target_allele"
        return result
    try:
        remaining_time(deadline)
        if not identity_anchor:
            reference = fasta.fetch(canonical, target["pos"] - 1,
                                    target["pos"] - 1 + len(target["ref"])).upper()
            if reference != target["ref"]:
                raise AlleleResolutionUnavailable("target_reference_sequence_mismatch")
        aliases = {row[0] for row in connection.execute(
            "select alias from chromosome_aliases where canonical=?", (canonical,)
        )}
        rows = {}
        source_loci = set()
        # Normalized coordinates can precede original source coordinates in a repeat.
        # Fetch all source loci that establish the exact identity, plus overlapping calls
        # at each original locus. A shifted exact allele does not require a separate
        # reference assertion at its normalized coordinate. Failed duplicates stay visible.
        positions = set() if hits else {target["pos"]}
        for hit in hits:
            remaining_time(deadline)
            row = connection.execute("select rowid as record_rowid,* from records where rowid=?",
                                     (hit["record_rowid"],)).fetchone()
            if row is None:
                raise AlleleResolutionUnavailable("normalized_source_record_missing")
            indexed = dict(row)
            source_loci.add((canonical, indexed["pos"], str(indexed["ref"]).upper()))
            positions.add(indexed["pos"])
            # Inspect starts within the full original REF span, even when repeat
            # normalization moved the public coordinate left of that span.
            interior = connection.execute("""select distinct r.pos from chromosome_aliases ca
                cross join records r on r.chrom=ca.alias where ca.canonical=?
                and r.pos>? and r.pos<=? limit ?""", (canonical, indexed["pos"],
                indexed["pos"] + len(indexed["ref"]) - 1, MAX_SOURCE_RECORDS+1)).fetchall()
            if len(interior) > MAX_SOURCE_RECORDS:
                raise AlleleResolutionUnavailable("source_call_limit")
            positions.update(row[0] for row in interior)
            if indexed["pos"] <= target["pos"] < indexed["pos"] + len(indexed["ref"]):
                positions.add(target["pos"])
        for position in sorted(positions):
            remaining_time(deadline)
            key = canonical, position
            if key not in cache:
                cache[key] = _records_at(
                    connection, canonical, position, array=array, source=source,
                    aliases=aliases, deadline=deadline,
                )
            for record in cache[key]:
                rows[(record["offset"], record.get("sample_index", 0))] = record
                if len(rows) > MAX_SOURCE_RECORDS:
                    raise AlleleResolutionUnavailable("source_call_limit")
        records = list(rows.values())
        result["source_records"] = records
        if not array:
            for record in records:
                remaining_time(deadline)
                anchor = fasta.fetch(canonical, record["pos"]-1,
                                     record["pos"]-1+len(record["ref"])).upper()
                if anchor != record["ref"].upper():
                    raise AlleleResolutionUnavailable("source_reference_sequence_mismatch")
        assessment = _combined_assessment(
            records, array=array, probe_id=target.get("rsid"), normalized_source_loci=source_loci,
        )
        result["call_assessment"] = assessment
        if not records:
            result.update(status="not_observed" if source_format == "vcf" else "not_covered",
                          reference_reason="plain_vcf_absence_is_inconclusive" if source_format == "vcf"
                          else "position_not_represented_in_complete_index")
            return result
        if assessment["status"] != "called":
            result.update(
                status="no_call" if assessment["status"] == "no_call" else "unavailable",
                reference_reason=assessment.get("reason") or "source_genotype_" + assessment["status"],
            )
            return result
        if identity_anchor:
            # A source identifier locates its original call without resolving a public
            # allele. REF==ALT is an internal anchor, never a reference/absence claim.
            result.update(status="not_observed", exact_match=False, target_dosage=None,
                          reference_reason="public_allele_mapping_unresolved")
            return result
        if hits:
            dosage = _target_dosage(hits, records, assessment, target, array=array)
            result["target_dosage"] = dosage
            result.update(status="observed" if dosage is not None and dosage > 0 else "unavailable",
                          reference_reason=None if dosage is not None and dosage > 0 else "target_dosage_unresolved")
            return result
        if array:
            # Zero inference requires the requested biallelic model. Literal
            # positive array observations, including third alleles, are handled above.
            if (
                len(target["ref"]) != 1
                or len(target["alt"]) != 1
                or target["ref"] not in "ACGT"
                or target["alt"] not in "ACGT"
                or any(base not in {target["ref"], target["alt"]} for base in assessment["alleles"])
            ):
                result["reference_reason"] = "array_allele_model_unsupported"
            else:
                assayed_reference = all(base == target["ref"] for base in assessment["alleles"])
                result.update(status="reference_confirmed" if assayed_reference else "not_observed",
                              target_dosage=0, reference_scope="assayed_array_marker",
                              reference_reason="array_sequencing_metrics_unavailable")
            return result
        if all(not record.get("is_variant") for record in records):
            status, reason = _reference_status(
                records, assessment, target, canonical=canonical, source_format=source_format, fasta=fasta,
            )
            result.update(status=status, reference_reason=reason,
                          target_dosage=0 if status == "reference_confirmed" else None)
            return result
        # A complete simple source genotype can report another ALT at this site.
        # Projection of an MNV or indel cannot supply target dosage or absence.
        if len(target["ref"]) == len(target["alt"]) == 1 and all(
               record["pos"] == target["pos"] and record["ref"].upper() == target["ref"]
               and len(record["ref"]) == 1
               and all(len(allele) == 1 for allele in assessment["alleles"])
               for record in records):
            if assessment.get("quality") == "pass":
                result.update(status="not_observed", target_dosage=0)
            else:
                reasons = assessment.get("quality_reasons", [])
                result["reference_reason"] = reasons[0] if reasons else "source_quality_insufficient"
        else:
            result["reference_reason"] = "sequence_allele_not_identical"
            projection = _site_projection(records, target)
            if projection is not None:
                result["site_projection"] = projection
        return result
    except AlleleResolutionUnavailable as error:
        result["reference_reason"] = error.reason
        return result
    except TimeoutError:
        raise
    except (OSError, UnicodeError, ValueError, KeyError, TypeError):
        result["reference_reason"] = "source_record_unreadable"
        return result


def resolve_allele_targets(reader: ActiveGenomeIndexReader, targets: list[dict], *,
                           reference_fasta: Path, expected_reference_sha256: str,
                           deadline: float, observed_only: bool = False) -> dict:
    """Resolve target order with one indexed exact join and bounded original assessment.

    observed_only accepts large public gene inventories but resolves only exact observed
    candidates. It does not perform exhaustive reference checks across a gene. Normal
    mode resolves every mapped rsID allele, including missing sites and reference calls.
    The expected digest comes from a trusted reference receipt; it must equal the
    index metadata. No clinical snapshot is read here, and no source/index is modified.
    relevant_target_ordinals is internal, complete exact identity inventory metadata.
    """
    from pysam import FastaFile

    if (not isinstance(targets, list) or len(targets) > MAX_INPUT_TARGETS
            or not math.isfinite(deadline) or not isinstance(observed_only, bool)
            or not isinstance(expected_reference_sha256, str)
            or not re.fullmatch(r"[0-9a-f]{64}", expected_reference_sha256)):
        raise ValueError("invalid_allele_resolution_request")
    prepared = []
    for target in targets:
        remaining_time(deadline)
        if (not isinstance(target, dict) or not all(isinstance(target.get(field), str) and target[field]
                                                   for field in ("chrom", "ref", "alt"))
                or type(target.get("pos")) is not int or target["pos"] <= 0
                or ("rsid" in target and target["rsid"] is not None and not isinstance(target["rsid"], str))):
            raise ValueError("invalid_allele_resolution_target")
        prepared.append({**target, "ref": target["ref"].upper(), "alt": target["alt"].upper()})
    reader.ensure_ready()
    if not reader.complete:
        raise AlleleResolutionUnavailable("active_genome_index_incomplete")
    with ExitStack() as stack:
        connection = stack.enter_context(reader.connect())
        from .normalized_alleles import normalized_allele_index_metadata
        metadata = _metadata(connection)
        normalized_allele_index_metadata(connection)
        _index_identity(metadata, reader.genome_build or metadata.get("genome_build"),
                        expected_reference_sha256)
        source_format = metadata.get("source_format")
        if source_format not in ARRAY_FORMATS | {"vcf", "gvcf"}:
            raise AlleleResolutionUnavailable("unsupported_source_assessment")
        # Main is read-only. TEMP writes do not change the managed genome index.
        connection.set_progress_handler(lambda: int(time.monotonic() >= deadline), 1000)
        connection.execute("create temp table resolution_targets(target_id integer primary key,chrom text,pos integer,ref text,alt text)")
        for start in range(0, len(prepared), 1000):
            remaining_time(deadline)
            connection.executemany("insert into temp.resolution_targets values(?,?,?,?,?)", [
                (ordinal, target["chrom"], target["pos"], target["ref"], target["alt"])
                for ordinal, target in enumerate(prepared[start:start + 1000], start=start)
            ])
        try:
            # Identify the complete inventory before the source-assessment bound.
            # EXISTS avoids materializing every duplicate original row for 100k targets.
            by_target = set()
            for kind, index, ref_clause in (
                ("variant_call", "normalized_seq_idx", "and k.ref=t.ref"),
                ("array_call", "normalized_array_idx", "and length(t.ref)=1 and t.ref in ('A','C','G','T') and length(t.alt)=1 and t.alt in ('A','C','G','T')"),
            ):
                sql = f"""select distinct t.target_id from temp.resolution_targets t
                    cross join chromosome_aliases ca on ca.alias=t.chrom
                    where exists(select 1 from normalized_alleles k indexed by {index}
                        where k.kind='{kind}' and k.chrom=ca.canonical
                          and k.pos=t.pos and k.alt=t.alt {ref_clause})"""
                for row in connection.execute(sql):
                    remaining_time(deadline)
                    by_target.add(row[0])
            relevant = sorted(by_target) if observed_only else list(range(len(prepared)))
            reference_fasta = Path(reference_fasta)
            index = Path(str(reference_fasta) + ".fai")
            if (reference_fasta.is_symlink() or index.is_symlink()
                    or not reference_fasta.is_file() or not index.is_file()):
                raise AlleleResolutionUnavailable("reference_context_unavailable")
            fasta = stack.enter_context(FastaFile(str(reference_fasta), filepath_index=str(index)))
            source = None if source_format in ARRAY_FORMATS else _OriginalSource(reader, metadata, stack)
            aliases = dict(connection.execute("select alias,canonical from chromosome_aliases"))
            cache = {}
            results = []
            unavailable_reason = "resolution_target_limit" if len(relevant) > MAX_RESOLVED_TARGETS else None
            try:
                for ordinal in relevant[:MAX_RESOLVED_TARGETS]:
                    remaining_time(deadline)
                    target = prepared[ordinal]
                    kind = "array_call" if source_format in ARRAY_FORMATS else "variant_call"
                    index = "normalized_array_idx" if kind == "array_call" else "normalized_seq_idx"
                    ref_clause = "" if kind == "array_call" else "and k.ref=?"
                    parameters = [target["chrom"], target["pos"], target["alt"]]
                    if ref_clause:
                        parameters.append(target["ref"])
                    hits = [dict(row) for row in connection.execute(f"""select k.* from chromosome_aliases ca
                        cross join normalized_alleles k indexed by {index}
                        where ca.alias=? and k.chrom=ca.canonical and k.kind='{kind}'
                        and k.pos=? and k.alt=? {ref_clause} limit {MAX_SOURCE_RECORDS+1}""", parameters)] if ordinal in by_target else []
                    result = _resolve_target(
                        connection, target, hits,
                        canonical=aliases.get(prepared[ordinal]["chrom"]), source_format=source_format,
                        source=source, fasta=fasta, deadline=deadline, cache=cache,
                    )
                    results.append({**result, "target_ordinal": ordinal})
            except (TimeoutError, sqlite3.OperationalError):
                if time.monotonic() < deadline:
                    raise
                unavailable_reason = "allele_resolution_deadline_exceeded"
            omitted = len(relevant) - len(results)
            return {
                "results": results, "total_target_count": len(prepared),
                "relevant_target_count": len(relevant), "returned_target_count": len(results),
                "omitted_target_count": omitted, "relevant_target_ordinals": sorted(by_target),
                "status": "partial" if omitted else "complete",
                **({"unavailable_reason": unavailable_reason} if unavailable_reason else {}),
            }
        except sqlite3.OperationalError:
            remaining_time(deadline)
            raise AlleleResolutionUnavailable("allele_resolution_database_unavailable") from None
        finally:
            connection.set_progress_handler(None, 0)
