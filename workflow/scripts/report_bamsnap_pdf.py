import glob
import logging
import os
import re

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)

from PIL import Image, ImageDraw, ImageFont

bamsnap_dir = snakemake.input.bamsnap_dir
output_pdf = snakemake.output.pdf


def _sort_key(path):
    """Sort PNG files by chromosome (numeric) then position."""
    name = os.path.splitext(os.path.basename(path))[0]  # e.g. "chr7_55174772"
    parts = name.rsplit("_", 1)
    chrom = parts[0]
    pos = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
    suffix = chrom.lstrip("chr").lstrip("Chr")
    try:
        return (int(suffix), pos)
    except ValueError:
        order = {"X": 100, "Y": 101, "M": 102, "MT": 102}
        return (order.get(suffix, 200), pos)


images_dir = os.path.join(bamsnap_dir, "bamsnap_images")
png_files = sorted(glob.glob(os.path.join(images_dir, "*.png")), key=_sort_key)

if not png_files:
    logging.warning(f"No PNG files found in {images_dir}; writing placeholder PDF")
    img = Image.new("RGB", (800, 120), "white")
    draw = ImageDraw.Draw(img)
    draw.text((20, 50), "No bamsnap images available.", fill="black")
    img.save(output_pdf, "PDF")
else:
    pages = [Image.open(f).convert("RGB") for f in png_files]
    logging.info(f"Writing {len(pages)} images to {output_pdf}")
    pages[0].save(output_pdf, "PDF", save_all=True, append_images=pages[1:])
    logging.info("Done")
