from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import Any

from .record_kinds import RECORD_KIND_ARRAY_CALL, RECORD_KIND_ARRAY_NO_CALL

JsonObject = dict[str, Any]

ARRAY_NO_CALLS = {"", ".", "--", "00", "NN"}
SUPPORTED_ARRAY_BASES = {"A", "C", "G", "T", "I", "D"}


@dataclass(frozen=True)
class ArrayGenotypeCall:
    raw_genotype: Any
    alleles: tuple[str, ...]
    status: str
    snp_eligible: bool


def classify_array_genotype(genotype: Any) -> ArrayGenotypeCall:
    """Keep the complete source call: a missing or complex token vetoes SNP evidence."""
    text = str(genotype or "").strip().upper()
    if not text or text in ARRAY_NO_CALLS:
        return ArrayGenotypeCall(genotype, (text,) if text else (), "no_call", False)
    if "/" in text or "|" in text:
        alleles = tuple(re.split(r"[/|]", text))
    elif len(text) in {1, 2} and all(base in SUPPORTED_ARRAY_BASES | {"N", "0", ".", "-"} for base in text):
        alleles = tuple(text)
    else:
        alleles = (text,)
    if any(allele in {"", ".", "-", "--", "0", "00", "N", "NN"} for allele in alleles):
        return ArrayGenotypeCall(genotype, alleles, "no_call", False)
    complete = 1 <= len(alleles) <= 2 and all(allele in SUPPORTED_ARRAY_BASES for allele in alleles)
    status = "called" if complete else "unresolved"
    snp_eligible = complete and all(allele in "ACGT" for allele in alleles)
    return ArrayGenotypeCall(genotype, alleles, status, snp_eligible)


def array_snp_call(genotype: Any, observed_alleles: Any) -> bool:
    call = classify_array_genotype(genotype)
    return call.snp_eligible and tuple(_observed_alleles(observed_alleles)) == call.alleles


def is_array_genotype_record(record: JsonObject) -> bool:
    return record.get("record_kind") in {RECORD_KIND_ARRAY_CALL, RECORD_KIND_ARRAY_NO_CALL}


def array_genotype_bases(record: JsonObject) -> list[str] | None:
    if record.get("record_kind") != RECORD_KIND_ARRAY_CALL:
        return None
    if str(record.get("filter") or "") not in {"PASS", "."}:
        return None
    call = classify_array_genotype(record.get("genotype"))
    if call.status != "called":
        return None
    bases = _observed_alleles(record.get("observed_alleles"))
    if not bases:
        return None
    if tuple(bases) != call.alleles:
        return None
    return bases


def called_genotype_tokens(genotype: Any) -> list[str]:
    text = str(genotype or "").strip().upper()
    if not text:
        return []
    if re.fullmatch(r"\d+(?:[/|]\d+)*", text):
        return re.split(r"[/|]", text)
    call = classify_array_genotype(genotype)
    return list(call.alleles) if call.status == "called" else []


def count_array_allele(
    record: JsonObject,
    *,
    target_allele: str,
    allowed_alleles: list[str] | tuple[str, ...] | set[str] | None = None,
) -> JsonObject:
    target = _single_array_allele(target_allele)
    if target is None:
        return {"status": "missing", "reason": "array_target_allele_not_single_base"}
    allowed: set[str] | None = None
    if allowed_alleles is not None:
        allowed = set()
        for allele in allowed_alleles:
            normalized = _single_array_allele(allele)
            if normalized is None:
                return {"status": "missing", "reason": "array_allele_model_not_single_base"}
            allowed.add(normalized)
        allowed.add(target)

    bases = array_genotype_bases(record)
    if bases is None:
        return {"status": "missing", "reason": "missing_genotype"}
    if allowed is not None and any(base not in allowed for base in bases):
        return {
            "status": "missing",
            "reason": "array_genotype_allele_outside_allowed_alleles",
            "allele_bases": bases,
        }
    return {
        "status": "matched",
        "allele_bases": bases,
        "dosage": float(sum(1 for base in bases if base == target)),
        "ploidy": len(bases),
    }


def _single_array_allele(value: str | None) -> str | None:
    allele = str(value or "").strip().upper()
    if len(allele) != 1 or allele not in SUPPORTED_ARRAY_BASES:
        return None
    return allele


def _observed_alleles(value: Any) -> list[str]:
    raw = value
    if isinstance(raw, str):
        try:
            raw = json.loads(raw)
        except json.JSONDecodeError:
            return []
    if not isinstance(raw, (list, tuple)):
        return []
    return [str(base).strip().upper() for base in raw if str(base).strip()]
