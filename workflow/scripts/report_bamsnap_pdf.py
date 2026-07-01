import glob
import logging
import os
import sys
import traceback

log_file = open(snakemake.log[0], "w")
sys.stderr = log_file

logging.basicConfig(
    stream=log_file,
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)

logging.info("report_bamsnap_pdf starting")
logging.info(f"input.bamsnap_dir = {snakemake.input.bamsnap_dir}")
logging.info(f"output.pdf        = {snakemake.output.pdf}")

try:
    from PIL import Image, ImageDraw
    logging.info(f"Pillow version: {Image.__version__}")
except ImportError as e:
    logging.error(f"Pillow not available: {e}")
    sys.exit(1)


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
    bamsnap_dir = snakemake.input.bamsnap_dir
    output_pdf = snakemake.output.pdf

    images_dir = os.path.join(bamsnap_dir, "bamsnap_images")
    logging.info(f"Looking for PNGs in {images_dir}")
    logging.info(f"images_dir exists: {os.path.isdir(images_dir)}")

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
        logging.info(f"Output created: {output_pdf} ({os.path.getsize(output_pdf)} bytes)")
    else:
        logging.error(f"Output file was NOT created: {output_pdf}")
        sys.exit(1)

except Exception:
    logging.error(traceback.format_exc())
    sys.exit(1)
finally:
    log_file.flush()
