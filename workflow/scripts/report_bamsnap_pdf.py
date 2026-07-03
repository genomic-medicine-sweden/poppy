#!/usr/bin/env python3
import glob
import logging
import os
import sys
import traceback

from PIL import Image, ImageDraw

bamsnap_dir = sys.argv[1]
output_pdf = sys.argv[2]

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)

logging.info(f"bamsnap_dir = {bamsnap_dir}")
logging.info(f"output_pdf = {output_pdf}")


def _sort_key(path):
    name = os.path.splitext(os.path.basename(path))[0]
    parts = name.rsplit("_", 1)
    chrom = parts[0]
    pos = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
    suffix = chrom.lstrip("chr").lstrip("Chr")
    try:
        return (int(suffix), pos)
    except ValueError:
        order = {"X": 100, "Y": 101, "M": 102, "MT": 102}
        return (order.get(suffix, 200), pos)


try:
    images_dir = os.path.join(bamsnap_dir, "bamsnap_images")
    logging.info(f"Looking for PNGs in {images_dir} (exists: {os.path.isdir(images_dir)})")

    png_files = sorted(glob.glob(os.path.join(images_dir, "*.png")), key=_sort_key)
    logging.info(f"Found {len(png_files)} PNG files")

    if not png_files:
        logging.warning("No PNG files found; writing placeholder PDF")
        img = Image.new("RGB", (800, 120), "white")
        draw = ImageDraw.Draw(img)
        draw.text((20, 50), "No bamsnap images available.")
        img.save(output_pdf, "PDF")
    else:
        pages = [Image.open(f).convert("RGB") for f in png_files]
        logging.info(f"Writing {len(pages)} pages to {output_pdf}")
        pages[0].save(output_pdf, "PDF", save_all=True, append_images=pages[1:])

    if os.path.isfile(output_pdf):
        logging.info(f"Done — {os.path.getsize(output_pdf)} bytes")
    else:
        logging.error(f"Output file not created: {output_pdf}")
        sys.exit(1)

except Exception:
    logging.error(traceback.format_exc())
    sys.exit(1)
