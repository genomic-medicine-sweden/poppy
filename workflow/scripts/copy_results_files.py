import subprocess


def main():
    with open(snakemake.log[0], "a") as logfile:
        for infile, outfile in zip(snakemake.input, snakemake.output):
            subprocess.run(
                ["rsync", "--update", "-a", infile, outfile],
                stdout=logfile,
                stderr=logfile,
            )


if __name__ == "__main__":
    main()
