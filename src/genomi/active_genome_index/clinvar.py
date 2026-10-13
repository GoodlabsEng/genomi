from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

from ..clinvar_match_model import (
    MATCH_BASIS_CONSUMER_ARRAY_ALLELE_INFERENCE,
    MATCH_BASIS_EXACT_ALLELE,
    MATCH_BASIS_MULTIALLELIC_ALT,
)
from .record_kinds import RECORD_KIND_ARRAY_CALL, RECORD_KIND_VARIANT_CALL
from .array_genotypes import array_snp_call
from .reader import ActiveGenomeIndexReader


JsonObject = dict[str, Any]

_CLINVAR_SELECTED_RECORD_COLUMN_NAMES = (
    "record_rowid",
    "chrom",
    "chrom_sort",
    "pos",
    "rsid",
    "ref",
    "alt",
    "qual",
    "filter",
    "info",
    "sample_index",
    "sample_name",
    "format",
    "genotype",
    "depth",
    "genotype_quality",
    "record_kind",
    "observed_alleles",
    "clinvar_batch_id",
    "clinvar_match_mode",
    "clinvar_match_alt",
)
_CLINVAR_AGI_RECORD_COLUMN_NAMES = _CLINVAR_SELECTED_RECORD_COLUMN_NAMES[:-3]
_CLINVAR_SELECTED_RECORD_COLUMNS = ", ".join(_CLINVAR_SELECTED_RECORD_COLUMN_NAMES)
_CLINVAR_AGI_RECORD_COLUMNS = ", ".join(_CLINVAR_AGI_RECORD_COLUMN_NAMES)


def selected_clinvar_record_columns(alias: str | None = None) -> str:
    if alias is None:
        return _CLINVAR_SELECTED_RECORD_COLUMNS
    return ", ".join(f"{alias}.{column}" for column in _CLINVAR_SELECTED_RECORD_COLUMN_NAMES)


def stage_clinvar_match_records(
    reader: ActiveGenomeIndexReader,
    connection: sqlite3.Connection,
) -> JsonObject:
    """Materialize the AGI rows ClinVar matching may read into temp tables."""

    reader.ensure_ready()
    connection.create_function("array_snp_call", 2, array_snp_call, deterministic=True)
    alias = "_agi_clinvar_records"
    reader.attach_to(connection, alias)
    try:
        _ensure_ready_for_clinvar_match(connection, alias, reader.agi_path)
        connection.executescript(
            """
            drop table if exists temp.selected_active_genome_index_agi_records;
            drop table if exists temp.selected_active_genome_index_records;
            create temp table selected_active_genome_index_agi_records (
                record_rowid integer not null,
                chrom text not null,
                chrom_sort integer,
                pos integer not null,
                rsid text,
                ref text,
                alt text,
                qual text,
                filter text,
                info text,
                sample_index integer,
                sample_name text,
                format text,
                genotype text,
                depth integer,
                genotype_quality integer,
                record_kind text,
                observed_alleles text
            );
            create temp table selected_active_genome_index_records (
                record_rowid integer not null,
                chrom text not null,
                chrom_sort integer,
                pos integer not null,
                rsid text,
                ref text,
                alt text,
                qual text,
                filter text,
                info text,
                sample_index integer,
                sample_name text,
                format text,
                genotype text,
                depth integer,
                genotype_quality integer,
                record_kind text,
                observed_alleles text,
                clinvar_batch_id text not null,
                clinvar_match_mode text not null,
                clinvar_match_alt text not null
            );
            create index selected_active_genome_index_agi_records_order_idx
                on selected_active_genome_index_agi_records(chrom_sort, pos, record_rowid, sample_index);
            create index selected_active_genome_index_records_locus_idx
                on selected_active_genome_index_records(chrom, pos);
            """
        )
        sql = f"""
            insert into temp.selected_active_genome_index_agi_records (
                record_rowid, chrom, chrom_sort, pos, rsid, ref, alt, qual, filter,
                info, sample_index, sample_name, format, genotype, depth,
                genotype_quality, record_kind, observed_alleles
            )
            select rowid, chrom, chrom_sort, pos, rsid, ref, alt, qual, filter,
                   info, sample_index, sample_name, format, genotype, depth,
                   genotype_quality, record_kind, observed_alleles
            from {alias}.records
            where record_kind in ('{RECORD_KIND_VARIANT_CALL}', '{RECORD_KIND_ARRAY_CALL}')
        """
        connection.execute(sql)
        connection.execute(
            """
            insert into temp.selected_active_genome_index_records (
                record_rowid, chrom, chrom_sort, pos, rsid, ref, alt, qual, filter,
                info, sample_index, sample_name, format, genotype, depth,
                genotype_quality, record_kind, observed_alleles,
                clinvar_batch_id, clinvar_match_mode, clinvar_match_alt
            )
            select distinct
                   r.record_rowid, r.chrom, r.chrom_sort, r.pos, r.rsid, r.ref, r.alt,
                   r.qual, r.filter, r.info, r.sample_index, r.sample_name, r.format,
                   r.genotype, r.depth, r.genotype_quality, r.record_kind, r.observed_alleles,
                   case
                       when instr(r.alt, ',') > 0
                           then cast(r.record_rowid as text) || ':' || upper(observed.value)
                       else cast(r.record_rowid as text)
                   end as clinvar_batch_id,
                   case
                       when instr(r.alt, ',') > 0 then ?
                       else ?
                   end as clinvar_match_mode,
                   upper(observed.value) as clinvar_match_alt
            from temp.selected_active_genome_index_agi_records as r
            join json_each(r.observed_alleles) as observed
            where r.record_kind = ?
              and r.alt not in ('', '.')
              and upper(observed.value) <> upper(r.ref)
              and instr(',' || upper(r.alt) || ',', ',' || upper(observed.value) || ',') > 0
            """,
            (
                MATCH_BASIS_MULTIALLELIC_ALT,
                MATCH_BASIS_EXACT_ALLELE,
                RECORD_KIND_VARIANT_CALL,
            ),
        )
        connection.execute(
            """
            insert into temp.selected_active_genome_index_records (
                record_rowid, chrom, chrom_sort, pos, rsid, ref, alt, qual, filter,
                info, sample_index, sample_name, format, genotype, depth,
                genotype_quality, record_kind, observed_alleles,
                clinvar_batch_id, clinvar_match_mode, clinvar_match_alt
            )
            select distinct
                   r.record_rowid, r.chrom, r.chrom_sort, r.pos, r.rsid, r.ref, r.alt,
                   r.qual, r.filter, r.info, r.sample_index, r.sample_name, r.format,
                   r.genotype, r.depth, r.genotype_quality, r.record_kind, r.observed_alleles,
                   cast(r.record_rowid as text) || ':' || upper(observed.value) || ':array'
                       as clinvar_batch_id,
                   ? as clinvar_match_mode,
                   upper(observed.value) as clinvar_match_alt
            from temp.selected_active_genome_index_agi_records as r
            join json_each(r.observed_alleles) as observed
            where r.record_kind = ?
              and array_snp_call(r.genotype,r.observed_alleles)
              and upper(observed.value) in ('A', 'C', 'G', 'T')
            """,
            (
                MATCH_BASIS_CONSUMER_ARRAY_ALLELE_INFERENCE,
                RECORD_KIND_ARRAY_CALL,
            ),
        )
        connection.commit()
        return {"source_format": _source_format(connection, alias)}
    finally:
        if connection.in_transaction:
            connection.rollback()
        connection.execute(f"detach database {alias}")


def selected_clinvar_match_records_cte_sql(
    *,
    pass_only: bool,
    max_records: int | None,
    cross_build: bool = False,
) -> str:
    if cross_build:
        return """
            with selected_records as (
                select """ + _CLINVAR_SELECTED_RECORD_COLUMNS + """,
                       sample_chrom_original, sample_pos_original
                from temp.lifted_selected_records
            )
        """
    agi_cte = _clinvar_agi_records_cte_sql(max_records=max_records)
    sql = """
            with """ + agi_cte + """,
            selected_agi_records as (
                select *
                from agi_records
                where 1 = 1
        """
    if pass_only:
        sql += " and filter in ('PASS', '.')"
    sql += """
            ),
            selected_records as (
                select """ + selected_clinvar_record_columns("a") + """
                from temp.selected_active_genome_index_records as a
                join selected_agi_records as s
                  on s.record_rowid = a.record_rowid
            )
        """
    return sql


def count_selected_clinvar_non_pass_records(
    connection: sqlite3.Connection,
    *,
    pass_only: bool,
    max_records: int | None,
) -> int:
    if not pass_only:
        return 0
    selection_params = (max_records,) if max_records is not None else ()
    row = connection.execute(
        f"""
        with {_clinvar_agi_records_cte_sql(max_records=max_records)}
        select count(*) as skipped_non_pass
        from agi_records
        where filter not in ('PASS', '.')
        """,
        selection_params,
    ).fetchone()
    return int(row["skipped_non_pass"] or 0)


def selected_clinvar_chrom_style(
    connection: sqlite3.Connection,
    *,
    pass_only: bool,
    max_records: int | None,
    cross_build: bool = False,
) -> str:
    selection_params: tuple[Any, ...] = (
        () if cross_build else ((max_records,) if max_records is not None else ())
    )
    row = connection.execute(
        f"""
        {selected_clinvar_match_records_cte_sql(pass_only=pass_only, max_records=max_records, cross_build=cross_build)}
        select
            coalesce(sum(case when chrom like 'chr%' then 1 else 0 end), 0) as chr_rows,
            count(*) as total_rows
        from selected_records
        """,
        selection_params,
    ).fetchone()
    return _chrom_style_from_counts(int(row["chr_rows"]), int(row["total_rows"]))


def populate_lifted_clinvar_match_records_table(
    connection: sqlite3.Connection,
    *,
    source_build: str,
    target_build: str,
    pass_only: bool,
    max_records: int | None,
) -> tuple[int, int]:
    from ..runtime.liftover import get_liftover

    lifter = get_liftover(source_build, target_build)
    connection.executescript(
        """
        drop table if exists temp.lifted_selected_records;
        create temp table lifted_selected_records (
            record_rowid integer not null,
            sample_chrom_original text not null,
            sample_pos_original integer not null,
            chrom text not null,
            chrom_sort integer,
            pos integer not null,
            rsid text,
            ref text,
            alt text,
            info text,
            qual text,
            filter text,
            sample_index integer,
            sample_name text,
            format text,
            genotype text,
            depth integer,
            genotype_quality integer,
            record_kind text,
            observed_alleles text,
            clinvar_batch_id text not null,
            clinvar_match_mode text not null,
            clinvar_match_alt text not null
        );
        create index lifted_selected_records_locus_idx
            on lifted_selected_records(chrom, pos);
        """
    )
    selection_params = (max_records,) if max_records is not None else ()
    source_rows = connection.execute(
        f"""
        {selected_clinvar_match_records_cte_sql(pass_only=pass_only, max_records=max_records)}
        select * from selected_records
        """,
        selection_params,
    ).fetchall()

    lifted_alleles = 0
    dropped_alleles = 0
    insert_buffer: list[tuple[Any, ...]] = []
    for row in source_rows:
        sample_chrom = row["chrom"]
        sample_pos = int(row["pos"])
        result = lifter.lift_position_full(sample_chrom, sample_pos)
        if result is None or result[2] != "+":
            dropped_alleles += 1
            continue
        lifted_chrom, lifted_pos, _strand = result
        insert_buffer.append(
            (
                int(row["record_rowid"]),
                sample_chrom,
                sample_pos,
                lifted_chrom,
                int(row["chrom_sort"]) if row["chrom_sort"] is not None else None,
                lifted_pos,
                row["rsid"],
                row["ref"],
                row["alt"],
                row["info"],
                row["qual"],
                row["filter"],
                row["sample_index"],
                row["sample_name"],
                row["format"],
                row["genotype"],
                row["depth"],
                row["genotype_quality"],
                row["record_kind"],
                row["observed_alleles"],
                row["clinvar_batch_id"],
                row["clinvar_match_mode"],
                row["clinvar_match_alt"],
            )
        )
        lifted_alleles += 1
    if insert_buffer:
        connection.executemany(
            """
            insert into temp.lifted_selected_records (
                record_rowid, sample_chrom_original, sample_pos_original,
                chrom, chrom_sort, pos, rsid, ref, alt, info, qual, filter,
                sample_index, sample_name, format, genotype, depth, genotype_quality, record_kind, observed_alleles,
                clinvar_batch_id, clinvar_match_mode, clinvar_match_alt
            )
            values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            insert_buffer,
        )
    return lifted_alleles, dropped_alleles


def _clinvar_agi_records_cte_sql(*, max_records: int | None) -> str:
    sql = """
            agi_records as (
                select """ + _CLINVAR_AGI_RECORD_COLUMNS + """
                from temp.selected_active_genome_index_agi_records
                where record_kind in ('variant_call', 'array_call')
        """
    if max_records is not None:
        sql += " order by chrom_sort, pos, record_rowid, sample_index"
        sql += " limit ?"
    sql += ")"
    return sql


def _chrom_style_from_counts(chr_rows: int, total_rows: int) -> str:
    if total_rows <= 0:
        return "empty"
    if chr_rows <= 0:
        return "bare"
    if chr_rows == total_rows:
        return "chr"
    return "mixed"


def _ensure_ready_for_clinvar_match(connection: sqlite3.Connection, alias: str, agi_path: Path) -> None:
    stats_count = connection.execute(f"select count(*) from {alias}.stats").fetchone()[0]
    index_names = {
        str(row["name"])
        for row in connection.execute(
            f"""
            select name
            from {alias}.sqlite_master
            where type = 'index' and tbl_name = 'records'
            """
        )
    }
    required_indexes = {"records_export_idx", "records_variant_idx"}
    missing_indexes = sorted(required_indexes - index_names)
    if stats_count == 0 or missing_indexes:
        details = []
        if stats_count == 0:
            details.append("missing stats rows")
        if missing_indexes:
            details.append(f"missing query indexes: {', '.join(missing_indexes)}")
        raise RuntimeError(
            f"Active Genome Index is incomplete for ClinVar refresh ({agi_path}): "
            f"{'; '.join(details)}. Rebuild the Active Genome Index from the source genome file once."
        )


def _source_format(connection: sqlite3.Connection, alias: str) -> str | None:
    try:
        row = connection.execute(f"select value from {alias}.metadata where key = 'source_format'").fetchone()
    except sqlite3.Error:
        return None
    if row is None:
        return None
    try:
        parsed = json.loads(str(row["value"]))
    except (TypeError, json.JSONDecodeError):
        parsed = row["value"]
    return str(parsed) if parsed else None
