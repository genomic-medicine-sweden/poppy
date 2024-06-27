
import gzip
import statistics
from pysam import VariantFile

vcf_files = snakemake.input.vcfs
artifact_panel = open(snakemake.output.artifact_panel, "w")

artefact_call_dict = {}
for file_name in vcf_files:
    for record in VariantFile(file_name).fetch():
        key = record.contig + "_" + str(record.pos) + "_" + record.info["SVTYPE"]
        af = float(record.info["AF"][0])

        if key not in artefact_call_dict:
            artefact_call_dict[key] = [af] # Ar forsta 0 i jonas kod antal? call_dict[caller][key]=[0,[]]
        else:
            artefact_call_dict[key].append(af)
            # vad gor FFPE_rm_dup_dict????


artifact_panel.write("Chromosome\tpos\tSV_type\tmedian_MAF\tsd_MAF\tnr_obs\n")
for key in artefact_call_dict:
    chrom = key.split("_")[0]
    pos = key.split("_")[1]
    svtype = key.split("_")[2]

    artefact_call_dict[key].sort() # blir de listan?
    median_af = statistics.median(artefact_call_dict[key])
    sd_af = 1000
    nr_obs = len(artefact_call_dict[key])
    if nr_obs >= 4:
        '''This is the sample variance s² with Bessel’s correction, also known as variance with N-1 degrees of freedom.
        Provided that the data points are representative (e.g. independent and identically distributed),
        the result should be an unbiased estimate of the true population variance.'''
        sd_af = statistics.stdev(artefact_call_dict[key])
    artifact_panel.write(chrom + "\t" + pos + "\t" + svtype + "\t" + str(median_af) + "\t" + str(sd_af) + "\t" + str(nr_obs) + "\n")