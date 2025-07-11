from collections import OrderedDict
import gzip
import logging
from pysam import VariantFile

from hydra_genetics.utils.io.chr import ChrTranslater
from hydra_genetics.utils.models.hotspot import MultiBpVariantData
from hydra_genetics.utils.models.hotspot import ReportClass
from hydra_genetics.utils.io.hotspot import Reader
from hydra_genetics.utils.io import utils

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def check_info_field(vcf_in, field="CSQ"):
    """
    Check whether the annotation field "CSQ" is in the VCF header.

    :param vcf_in: Input VCF file.
    :return: number of annotations in the CSQ field.
    """
    vcfobj = VariantFile(vcf_in, "r")
    header = vcfobj.header  # Requires to index the VCF file first if gzipped, set index file as input in smk rule
    if field in header.info.keys():
        logger.info(f"Field {field} is present in the VCF header.")
        return len(header.info[field].description.split("Format: ")[-1].split('|'))
    else:
        logger.warning(f"Field {field} is NOT present in the VCF header.")
        return 0


def add_missing_annotation(vcf_in, vcf_out, field="CSQ"):
    """
    Add missing annotation field to the VCF file and set all annotations to blank.

    :param vcf_in:
    :param vcf_out:
    :param field:
    :return:
    """
    nb_annot = check_info_field(vcf_in, field=field)
    if nb_annot:
        logger.info(f"Field {field} is present in the VCF header. All variants must have INFO/CSQ.")
        vcfobj = VariantFile(vcf_in, "r")
        variants = vcfobj.fetch()
        logger.info("Opening output vcf: {}".format(vcf_out))
        with VariantFile(vcf_out, 'w', header=vcfobj.header) as vcfobjout:
            for variant in variants:
                # print(f"{variant.chrom}\t{variant.pos}\t{variant.id}\t{variant.ref}\t{variant.alts}\t{variant.qual}")
                field_value = variant.info.get(field, None)
                if field_value is None:
                    logger.info(f"Field {field} is missing in variant {variant.chrom}:{variant.pos}. Adding it now.")
                    blank_csq = [""] * nb_annot
                    variant.info.update({"CSQ": "|".join(blank_csq)})
                    print(f"{variant.info[field]}")
                    print(variant.info.items())
                vcfobjout.write(variant)
        logger.info("Closing output vcf: {}".format(vcf_out))


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    annotated = check_info_field(snakemake.input.vcf, field=snakemake.params.field)
    if annotated:
        add_missing_annotation(snakemake.input.vcf, snakemake.output.vcf, snakemake.params.field)
