#!/usr/bin/env python3
from pathlib import Path
import logging


def find_metrics_header(lines):
    for idx, line in enumerate(lines):
        stripped = line.rstrip("\n")
        if stripped.startswith("MEDIAN_INSERT_SIZE\t") or stripped == "MEDIAN_INSERT_SIZE":
            return idx
    raise ValueError("Could not find Picard metrics header line starting with 'MEDIAN_INSERT_SIZE'.")


def find_first_data_row(lines, header_idx):
    for idx in range(header_idx + 1, len(lines)):
        stripped = lines[idx].strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue
        return idx
    raise ValueError("Could not find first data row after Picard metrics header.")


def normalize_insert_size(lines, min_insert_size, replace_insert_size):
    header_idx = find_metrics_header(lines)
    data_idx = find_first_data_row(lines, header_idx)
    header_line = lines[header_idx].rstrip("\n")
    data_line = lines[data_idx].rstrip("\n")
    headers = header_line.split("\t")
    values = data_line.split("\t")

    if "MEDIAN_INSERT_SIZE" not in headers:
        raise ValueError("Header line does not contain 'MEDIAN_INSERT_SIZE' column.")
    if "MEAN_INSERT_SIZE" not in headers:
        raise ValueError("Header line does not contain 'MEAN_INSERT_SIZE' column.")

    median_idx = headers.index("MEDIAN_INSERT_SIZE")
    mean_idx = headers.index("MEAN_INSERT_SIZE")

    if median_idx >= len(values):
        raise ValueError(
            "First data row after Picard metrics header does not contain enough columns "
            "to match the header."
        )

    original_value = values[median_idx]
    try:
        insert_size = float(original_value)
    except ValueError as exc:
        raise ValueError(
            f"MEDIAN_INSERT_SIZE value '{original_value}' is not numeric in first data row."
        ) from exc

    changed = False
    if insert_size < min_insert_size:
        values[median_idx] = str(replace_insert_size)
        values[mean_idx] = str(replace_insert_size)
        lines[data_idx] = "\t".join(values) + "\n"
        changed = True

    return lines, changed, original_value


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    min_insert_size = snakemake.params.min_insert_size
    replace_insert_size = snakemake.params.replace_insert_size

    input_path = Path(snakemake.input.metrics)
    output_path = Path(snakemake.output.metrics)

    logging.info(f"Reading insert size metrics from {input_path}")
    lines = input_path.read_text().splitlines(keepends=True)
    updated_lines, changed, original_value = normalize_insert_size(
        lines, min_insert_size, replace_insert_size
    )
    output_path.write_text("".join(updated_lines))

    if changed:
        logging.info(
            "MEDIAN_INSERT_SIZE was below %s (%s); adjusted to %s for pindel input.",
            min_insert_size,
            original_value,
            replace_insert_size,
        )
    else:
        logging.info("MEDIAN_INSERT_SIZE was %s; no adjustment needed.", original_value)
    logging.info(f"Wrote adjusted metrics to {output_path}")
