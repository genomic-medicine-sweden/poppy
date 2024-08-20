from pysam import VariantFile


def add_artifact_annotation_data(in_vcf_filename, artifacts_filename, out_vcf_filename):
    artifact_dict = {}
    with open(artifacts_filename, "r") as artifacts:
        next(artifacts)
        for lline in artifacts:
            line = lline.strip().split("\t")
            chrom = line[0]
            pos = line[1]
            type = line[2]
            median = line[3]
            sd = line[4]
            observations = line[5]

            artifact_dict[chrom + "_" + pos + "_" + type] = [median, sd, observations]

    # Create new vcf
    in_vcf = VariantFile(in_vcf_filename)
    new_header = in_vcf.header
    new_header.info.add("Artifact", "1", "String", "Number of observations of variant in panel samples")
    new_header.info.add("ArtifactMedian", "1", "String", "Artifact median MAFs in normal panel")
    new_header.info.add("ArtifactNrSD", "1", "String", "Number of Standard Deviations from artifacts in panel median")
    out_vcf = VariantFile(out_vcf_filename, "w", header=new_header)

    for record in in_vcf.fetch():
        key = record.contig + "_" + str(record.pos) + "_" + record.info["SVTYPE"]
        AF = float(record.info["AF"][0])

        Observations = 0
        Median = 0
        SD = 0
        if key in artifact_dict:
            Median = artifact_dict[key][0]
            SD = artifact_dict[key][1]
            Observations = artifact_dict[key][2]

        record.info["Artifact"] = str(Observations)
        record.info["ArtifactMedian"] = str(Median)

        nrsd = 1000
        if float(SD) == 1000.0 or float(SD) == 0.0:
            if float(Median) >= AF:
                nrsd = 0.0
            else:
                nrsd = 1000
        else:
            nrsd = (AF - float(Median)) / float(SD)

        record.info["ArtifactNrSD"] = str(nrsd)
        out_vcf.write(record)


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    add_artifact_annotation_data(
        snakemake.input.vcf,
        snakemake.input.artifacts,
        snakemake.output.vcf,
    )
