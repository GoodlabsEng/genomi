"""Bounded compatibility checks for headerless Nucleus-style SNP text, not vendor identity."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterator, TextIO

from ..array_genotypes import classify_array_genotype
from .text_io import _open_text_source

MARKER_DEFINITIONS = json.loads(Path(__file__).with_name("nucleus_text_markers.json").read_text())
_MARKERS = {marker["rsid"]: marker for marker in MARKER_DEFINITIONS["markers"]}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_CHROMOSOME = re.compile(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)")
_TOKEN = re.compile(r"(?:[ACGTNID]+|\.|-|--|0|00|NN|\*|<[A-Za-z0-9:_.*+\-]+>)", re.IGNORECASE)


def _invalid(reason: str) -> ValueError:
    return ValueError(f"nucleus_text_{reason}: Upload a VCF that declares its reference build and alleles.")


def iter_nucleus_text_rows(handle: TextIO) -> Iterator[dict]:
    """Every nonempty row must be four tab-separated fields and a complete token pair."""
    for index, raw in enumerate(handle):
        text = raw.rstrip("\r\n")
        if index == 0:
            text = text.removeprefix("\ufeff")
        if not text.strip():
            continue
        columns = text.split("\t")
        if len(columns) != 4:
            raise _invalid("layout_invalid")
        rsid, chrom, pos_text, genotype = columns
        tokens = genotype.split("/")
        if (not re.fullmatch(r"rs[1-9][0-9]*", rsid)
                or not _CHROMOSOME.fullmatch(chrom)
                or not re.fullmatch(r"[1-9][0-9]*", pos_text)
                or len(pos_text) > 16 or int(pos_text) > 2**53 - 1
                or len(tokens) != 2 or any(not _TOKEN.fullmatch(token) for token in tokens)):
            raise _invalid("layout_invalid")
        chrom = chrom.removeprefix("chr")
        if chrom == "M":
            chrom = "MT"
        yield {"rsid": rsid, "chrom": chrom, "pos": int(pos_text), "genotype": genotype}


def looks_like_nucleus_text(lines: list[str]) -> bool:
    first = next((line for line in lines if line.strip()), "").removeprefix("\ufeff")
    columns = first.split("\t")
    return len(columns) == 4 and columns[0].startswith("rs") and "/" in columns[3]


def validate_nucleus_text_source(source: Path, *, member_name: str | None = None) -> dict:
    """Public definitions stay local; all observed marker contradictions veto inference."""
    coordinate_markers: set[str] = set()
    coordinate_chromosomes: set[str] = set()
    orientation_markers: set[str] = set()
    orientation_chromosomes: set[str] = set()
    summary = {"total_rows": 0, "supported_single_base_calls": 0, "incomplete_calls": 0, "unsupported_calls": 0}
    with _open_text_source(source, member_name=member_name) as handle:
        for row in iter_nucleus_text_rows(handle):
            call = classify_array_genotype(row["genotype"])
            summary["total_rows"] += 1
            summary["supported_single_base_calls" if call.snp_eligible else
                    "incomplete_calls" if call.status == "no_call" else "unsupported_calls"] += 1
            marker = _MARKERS.get(row["rsid"])
            if marker is None:
                continue
            if (row["chrom"] != marker["chromosome"] or row["pos"] != marker["grch37Position"]
                    or marker["grch37Position"] == marker["grch38Position"]):
                raise _invalid("coordinate_conflict")
            coordinate_markers.add(row["rsid"])
            coordinate_chromosomes.add(row["chrom"])
            if not call.snp_eligible:
                continue
            forward = set(marker["forwardAlleles"])
            reverse = {allele.translate(_COMPLEMENT) for allele in forward}
            is_forward = set(call.alleles) <= forward
            is_reverse = set(call.alleles) <= reverse
            if not is_forward:
                raise _invalid("orientation_conflict" if is_reverse else "marker_alleles_unknown")
            if not is_reverse:
                orientation_markers.add(row["rsid"])
                orientation_chromosomes.add(row["chrom"])
    if (len(coordinate_markers) < MARKER_DEFINITIONS["minimumCoordinateMarkers"]
            or len(coordinate_chromosomes) < MARKER_DEFINITIONS["minimumCoordinateChromosomes"]
            or len(orientation_markers) < MARKER_DEFINITIONS["minimumOrientationMarkers"]
            or len(orientation_chromosomes) < MARKER_DEFINITIONS["minimumOrientationChromosomes"]):
        raise _invalid("evidence_insufficient")
    return {"reference_validation": {"markerSetVersion": MARKER_DEFINITIONS["version"],
                                    "coordinateMarkerCount": len(coordinate_markers),
                                    "orientationMarkerCount": len(orientation_markers)},
            "array_call_summary": summary}
