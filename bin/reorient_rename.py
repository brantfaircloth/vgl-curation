#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "biopython>=1.85",
#     "click>=8.1",
# ]
# ///

import sys

import click
from Bio import SeqIO


def read_config(path):
    """Parse a TSV of old_name, new_name, orientation (+ or -)."""
    entries = []
    seen = set()
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) == 2:
                old, new, orient = fields[0], fields[1], "+"
            elif len(fields) == 3:
                old, new, orient = fields
            else:
                raise click.ClickException(
                    f"{path}:{lineno}: expected 2 or 3 fields, got {len(fields)}"
                )
            if orient not in ("+", "-"):
                raise click.ClickException(
                    f"{path}:{lineno}: orientation must be + or -, got {orient!r}"
                )
            if old in seen:
                raise click.ClickException(
                    f"{path}:{lineno}: duplicate entry for {old}"
                )
            seen.add(old)
            entries.append((old, new, orient))
    return entries


@click.command()
@click.option(
    "-f",
    "--fasta",
    required=True,
    type=click.Path(exists=True),
    help="Input assembly FASTA",
)
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="TSV: old_name, new_name, orientation",
)
@click.option("-o", "--output", required=True, type=click.Path(), help="Output FASTA")
@click.option(
    "--drop-unlisted",
    is_flag=True,
    help="Omit sequences absent from the config (default: keep, in input order, after listed ones)",
)
@click.option(
    "--line-width", default=60, show_default=True, help="FASTA line wrap width"
)
def main(fasta, config, output, drop_unlisted, line_width):
    """Rename and reorient contigs according to a config file."""
    entries = read_config(config)
    index = SeqIO.index(fasta, "fasta")

    missing = [old for old, _, _ in entries if old not in index]
    if missing:
        raise click.ClickException(
            "config names absent from FASTA: " + ", ".join(missing)
        )

    listed = {old for old, _, _ in entries}
    n_rc = 0
    with open(output, "w") as out:
        writer = SeqIO.FastaIO.FastaWriter(out, wrap=line_width)
        for old, new, orient in entries:
            record = index[old]
            if orient == "-":
                record = record.reverse_complement()
                n_rc += 1
            record.id = new
            record.name = new
            record.description = ""
            writer.write_record(record)
        n_kept = 0
        if not drop_unlisted:
            for record in SeqIO.parse(fasta, "fasta"):
                if record.id in listed:
                    continue
                record.description = ""
                writer.write_record(record)
                n_kept += 1

    click.echo(
        f"{len(entries)} renamed ({n_rc} reverse-complemented), {n_kept} passed through",
        err=True,
    )


if __name__ == "__main__":
    sys.exit(main())
