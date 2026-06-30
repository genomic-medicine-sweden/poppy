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
            f"MEDIAN_INSERT_SIZE was below {min_insert_size} ({original_value}); "
            f"adjusted to {replace_insert_size} for pindel input."
        )
    else:
        logging.info(
            f"MEDIAN_INSERT_SIZE was {original_value}; no adjustment needed."
        )
    logging.info(f"Wrote adjusted metrics to {output_path}")
