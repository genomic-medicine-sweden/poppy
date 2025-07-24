# Changelog

## [0.2.1](https://github.com/genomic-medicine-sweden/poppy/compare/v0.2.0...v0.2.1) (2025-07-24)


### Bug Fixes

* write blank CSQ for pindel variants skipped by VEP ([#80](https://github.com/genomic-medicine-sweden/poppy/issues/80)) ([0b5be7a](https://github.com/genomic-medicine-sweden/poppy/commit/0b5be7a7d9843f4097a63e391c39b2a5ccebb1ed))

## [0.2.0](https://github.com/genomic-medicine-sweden/poppy/compare/v0.1.0...v0.2.0) (2024-12-13)


### Features

* add background to reference pipeline and annotation to snvs vcfs ([#70](https://github.com/genomic-medicine-sweden/poppy/issues/70)) ([f5c1903](https://github.com/genomic-medicine-sweden/poppy/commit/f5c1903686952f99439683850004391935ae6222))
* add pindel artifact annotation and filter ([#71](https://github.com/genomic-medicine-sweden/poppy/issues/71)) ([8aaa171](https://github.com/genomic-medicine-sweden/poppy/commit/8aaa171162433dd9f2f32fdd8266c1c2ba247005))
* update hydra genetics version to include software version into multiqc ([#63](https://github.com/genomic-medicine-sweden/poppy/issues/63)) ([19b33d3](https://github.com/genomic-medicine-sweden/poppy/commit/19b33d37d1bdf583f01d5663a61183a25b526cca))


### Bug Fixes

* add bam-files to purecn rules where bamlist is used ([#73](https://github.com/genomic-medicine-sweden/poppy/issues/73)) ([b657c84](https://github.com/genomic-medicine-sweden/poppy/commit/b657c849040ff1763771dfd1ab34e411dcd037f6))
* remove default string for all col when loading samples.tsv ([#67](https://github.com/genomic-medicine-sweden/poppy/issues/67)) ([9efae47](https://github.com/genomic-medicine-sweden/poppy/commit/9efae473e1cb1c9319e010f3ffe6aecd970cd4f8))
* remove pindel decompose and update filters ([#72](https://github.com/genomic-medicine-sweden/poppy/issues/72)) ([6efbfc0](https://github.com/genomic-medicine-sweden/poppy/commit/6efbfc088233e5947600d9331d4fd733c08bffdf))

## 0.1.0 (2024-03-08)


### Features

* add annotation and seq dict for pindel vcf ([046d2a5](https://github.com/genomic-medicine-sweden/poppy/commit/046d2a5c244d7b0d6689561e833ffc9fd5883f26))
* add artifafct annotation ([25c43a7](https://github.com/genomic-medicine-sweden/poppy/commit/25c43a7c4e8021f897b6b3f1499477e00ad5cf90))
* add CNV HTML report ([#40](https://github.com/genomic-medicine-sweden/poppy/issues/40)) ([373c012](https://github.com/genomic-medicine-sweden/poppy/commit/373c01292f80a6516b55f1e7ecb09c9623784500))
* add cnvkit rules ([2cce492](https://github.com/genomic-medicine-sweden/poppy/commit/2cce4921ba8fa0a8b2bb81803ab71e12d7888191))
* add cnvkit vcf to output ([fd40472](https://github.com/genomic-medicine-sweden/poppy/commit/fd40472493be30d5755872b1e027f12a011587a6))
* add GATK CNVs to SVDB output ([1698966](https://github.com/genomic-medicine-sweden/poppy/commit/1698966170331ea9ebf5492fabe91310586350f2))
* add GATK interval list to output files ([6829ca0](https://github.com/genomic-medicine-sweden/poppy/commit/6829ca07e1c08a6dbac8fe14137a0553f9eb9ccf))
* add GATK panel of normals to output ([c93d77e](https://github.com/genomic-medicine-sweden/poppy/commit/c93d77ec2df92509342fdde16e3223c9d7cce0ba))
* add lower ad filter for myd88 and cxcr4 (gbg) ([c0c8e62](https://github.com/genomic-medicine-sweden/poppy/commit/c0c8e6235505ec2e27c65093d39f6dc37760d8ea))
* add pindelfilter ([6b5d320](https://github.com/genomic-medicine-sweden/poppy/commit/6b5d320583916fa3cea253e85edf2b513205f511))
* add purecn ([#50](https://github.com/genomic-medicine-sweden/poppy/issues/50)) ([71bda05](https://github.com/genomic-medicine-sweden/poppy/commit/71bda0519d4aa65154dace50d0143f253d821fb0))
* add somatic soft filter vcf ([0dda471](https://github.com/genomic-medicine-sweden/poppy/commit/0dda471d4317c25fd18a8c2358eeb8945c866724))
* add svdb output ([11e9398](https://github.com/genomic-medicine-sweden/poppy/commit/11e939896b0551cc5c718a196521c69300104356))
* add vep annotation to pindel if not empty ([1d540ee](https://github.com/genomic-medicine-sweden/poppy/commit/1d540ee9175905694dd5ff4c4a05e11fbb1788c4))
* added pindel and individual vcf ([e0c7452](https://github.com/genomic-medicine-sweden/poppy/commit/e0c7452b0cb236f614c1b8d18c47663a0dba23d6))
* adding softwareversion files, and updating multiqc to latest to include said version log ([2cb422a](https://github.com/genomic-medicine-sweden/poppy/commit/2cb422ac8fc62b10747594e3c90ef698a0226b66))
* allow comments using "#" in samples and units ([6beeea8](https://github.com/genomic-medicine-sweden/poppy/commit/6beeea87cd249aec75a6022e12b1ead6b84646f6))
* annotated snv vcf ([7f4b3d9](https://github.com/genomic-medicine-sweden/poppy/commit/7f4b3d942d24c28bada10ab45eca6a4b0971ed18))
* better copy-number thresholds for hematology ([#53](https://github.com/genomic-medicine-sweden/poppy/issues/53)) ([c507302](https://github.com/genomic-medicine-sweden/poppy/commit/c50730206e0a4312762bff362e5f2b79fd3f770e))
* generate svdb database ([e0b545d](https://github.com/genomic-medicine-sweden/poppy/commit/e0b545d7a813e45de4a2fd5bb6923f9febb19ef1))
* integration test dataset ([186749a](https://github.com/genomic-medicine-sweden/poppy/commit/186749a7e8aba3741a8c21ff16a0558ac49dd001))
* make tumor content optional ([#54](https://github.com/genomic-medicine-sweden/poppy/issues/54)) ([b50a63c](https://github.com/genomic-medicine-sweden/poppy/commit/b50a63c28651a57114002f20a941b3114a3d67e7))
* more wildcard constraints ([58a3101](https://github.com/genomic-medicine-sweden/poppy/commit/58a3101e88256dc7868f366965551aa00af53f7c))
* normalize, decompose and add af and dp to pindel vcf ([38eda80](https://github.com/genomic-medicine-sweden/poppy/commit/38eda8068d0689385f66628254a9211ac2758910))
* reference pipeline improvements ([#43](https://github.com/genomic-medicine-sweden/poppy/issues/43)) ([407e1db](https://github.com/genomic-medicine-sweden/poppy/commit/407e1db96fa31c9baafba2eda02019b18f736702))
* specific freebayes filter for low qual ([0594e95](https://github.com/genomic-medicine-sweden/poppy/commit/0594e95212776732db65346ce9ad6efe1c1ea13c))
* specify pipeline output in yaml file ([59e94e9](https://github.com/genomic-medicine-sweden/poppy/commit/59e94e958688da4afc329465c43266c4a39f3d1e))
* supply config file on command line ([3c42289](https://github.com/genomic-medicine-sweden/poppy/commit/3c42289dd6c7f41a81d47a2e2b5785e0dd82d255))
* Update prealignment module to v0.1.0, as well as hydra-genetics to v0.9.1 ([8fa9003](https://github.com/genomic-medicine-sweden/poppy/commit/8fa9003cf00c5d9226cf72df00d3285bc7ba8c48))
* update reports module to v0.3.1 ([#45](https://github.com/genomic-medicine-sweden/poppy/issues/45)) ([43350b3](https://github.com/genomic-medicine-sweden/poppy/commit/43350b3830766f33ff501ac68e2223553368d218))
* update vepto 109 and hg38 ([2f65f7a](https://github.com/genomic-medicine-sweden/poppy/commit/2f65f7a18c869047d09d8e6b3e352c8971b5e16f))


### Bug Fixes

* add minimum hydra version ([59ae8ad](https://github.com/genomic-medicine-sweden/poppy/commit/59ae8ad319fb20f350e518c9c8e6cadcf7dc3a5f))
* add missing lines and chars in configs ([c89ea4c](https://github.com/genomic-medicine-sweden/poppy/commit/c89ea4c0c1c1bf815ed161470420c47d4c3bf462))
* add missing mapping bias to purecn and update cnvkit_batch inputname for normal pool ([8c4027e](https://github.com/genomic-medicine-sweden/poppy/commit/8c4027ed7030172246adf1741379f77fabcf062e))
* add softfilter flag ([428d4a2](https://github.com/genomic-medicine-sweden/poppy/commit/428d4a257baeb60e74673b632cb6d55840094aba))
* bug where files w/o wildcards were missing ([89f9e5f](https://github.com/genomic-medicine-sweden/poppy/commit/89f9e5fa36e52de25539cdcf60a8e2313acca87d))
* bump cnv_sv version ([87aae77](https://github.com/genomic-medicine-sweden/poppy/commit/87aae773cc0f7a278731e5d2fd4343b717eb4ec5))
* check that fastq files exist ([7120a2b](https://github.com/genomic-medicine-sweden/poppy/commit/7120a2bc65f7b7659ba7968444ecc6f1008c1a9b))
* **config:** correct fastqc path ([418d109](https://github.com/genomic-medicine-sweden/poppy/commit/418d1093cc1f7e08791f0109299472ef9e821b5a))
* containers and rules in config ([0accab4](https://github.com/genomic-medicine-sweden/poppy/commit/0accab48913f698a01e27c4db7cf61d63f897f25))
* ensemble vcf config ([a2b6ae1](https://github.com/genomic-medicine-sweden/poppy/commit/a2b6ae16c4fe98be5540921921d0e3f5f3114d74))
* fix container definition for `pindel_update_vcf` ([b50a63c](https://github.com/genomic-medicine-sweden/poppy/commit/b50a63c28651a57114002f20a941b3114a3d67e7))
* flesh out config validation schema ([9a76b67](https://github.com/genomic-medicine-sweden/poppy/commit/9a76b6794c1684a0846fe5426fa97b8a006e34a7))
* incorrect type for copy number thresholds ([#55](https://github.com/genomic-medicine-sweden/poppy/issues/55)) ([97cc0c7](https://github.com/genomic-medicine-sweden/poppy/commit/97cc0c76c63e336f9b2daccc94041e095c181a0b))
* indentation error in validation schema ([4ed18f3](https://github.com/genomic-medicine-sweden/poppy/commit/4ed18f33ff11cf39413747097e930a3f5a8711ba))
* make sure pindel if statement works on all files ([baeb542](https://github.com/genomic-medicine-sweden/poppy/commit/baeb542b6b8279f37adf9af05eafc12adc58ac3c))
* make the db_string config parameter mandatory ([39c314b](https://github.com/genomic-medicine-sweden/poppy/commit/39c314bbf5012ba05ea517969cdb21a51c488610))
* more concise validation error message ([1073900](https://github.com/genomic-medicine-sweden/poppy/commit/1073900b9216810cf342335cdd7af6589fda9302))
* pin pulp version for snakemake &lt;8.1.2 ([#47](https://github.com/genomic-medicine-sweden/poppy/issues/47)) ([c07beb7](https://github.com/genomic-medicine-sweden/poppy/commit/c07beb78567948f5a9852a19f1003882bbd2f37e))
* **reference:** allow for sample to be split over several lines in units_references.tsv ([e269995](https://github.com/genomic-medicine-sweden/poppy/commit/e269995db7d9f8e93063d765cd1a9c15afb5e1a6))
* set default reference name ([2981a9b](https://github.com/genomic-medicine-sweden/poppy/commit/2981a9bb414946cb312c6f9e832ae6ea16960061))
* spelling on softfilter filename ([9381cbb](https://github.com/genomic-medicine-sweden/poppy/commit/9381cbbd820e90e1750abbd6c0212ec8a64f8fd8))
* update minimum snakemake version ([4e90b8a](https://github.com/genomic-medicine-sweden/poppy/commit/4e90b8ae3412415e9177f9639adecc6dffe7efd5))
* update module versions ([eea7afe](https://github.com/genomic-medicine-sweden/poppy/commit/eea7afe8d1f12ac19e9422025f5b0e3bbda34d11))
* use repo config ([8492bd7](https://github.com/genomic-medicine-sweden/poppy/commit/8492bd77e87d54d52b1b2a9d6868b59ab0082051))
* **validation:** add cnvkit rules to schema ([f65b824](https://github.com/genomic-medicine-sweden/poppy/commit/f65b82485704b06adf4657d50726f821507acfee))


### Performance Improvements

* remove freebayes caller ([0217502](https://github.com/genomic-medicine-sweden/poppy/commit/02175027cb8a804d0b0015e9206f3905617330b5))


### Miscellaneous Chores

* release 0.1.0 ([29324f3](https://github.com/genomic-medicine-sweden/poppy/commit/29324f33f332b0e7a91e953ebec47017dc84e5c5))
* update to new develop for modules ([3c89f79](https://github.com/genomic-medicine-sweden/poppy/commit/3c89f79dcc2a6cb57a44b13f3b0d1e98adcf4493))
