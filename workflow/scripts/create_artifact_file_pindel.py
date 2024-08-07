import gzip
import statistics
from pysam import VariantFile

vcf_files = snakemake.input.vcfs
artifact_panel = open(snakemake.output.artifact_panel, "w")

artifact_call_dict = {}
for file_name in vcf_files:
    in_vcf = VariantFile(file_name)
    samplename = in_vcf.header.samples[0]
    for record in in_vcf.fetch():
        key = record.contig + "_" + str(record.pos) + "_" + record.info["SVTYPE"]
        af = float(record.info["AF"][0])

        if key not in artifact_call_dict:
            artifact_call_dict[key] = [[af], [samplename]]  # Ar forsta 0 i jonas kod antal? call_dict[caller][key]=[0,[]]
        else:
            artifact_call_dict[key][0].append(af)
            artifact_call_dict[key][1].append(samplename)
            # vad gor FFPE_rm_dup_dict????


artifact_panel.write("Chromosome\tpos\tSV_type\tmedian_MAF\tsd_MAF\tnr_obs\n")
for key in artifact_call_dict:
    chrom = key.split("_")[0]
    pos = key.split("_")[1]
    svtype = key.split("_")[2]

    artifact_call_dict[key][0].sort()
    artifact_call_dict[key][1] = list(dict.fromkeys(artifact_call_dict[key][1]))

    median_af = statistics.median(artifact_call_dict[key][0])
    sd_af = 1000
    nr_obs = len(artifact_call_dict[key][1])  # bara antal prover, varje prov kan ha flera.
    if nr_obs >= 4:
        """This is the sample variance s² with Bessel’s correction, also known as variance with N-1 degrees of freedom.
        Provided that the data points are representative (e.g. independent and identically distributed),
        the result should be an unbiased estimate of the true population variance."""
        sd_af = statistics.stdev(artifact_call_dict[key][0])
    artifact_panel.write(
        chrom + "\t" + pos + "\t" + svtype + "\t" + str(median_af) + "\t" + str(sd_af) + "\t" + str(nr_obs) + "\n"
    )
