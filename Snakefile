#!/usr/bin/env python3
import pathlib2

##needed to get BUSCO running in new folder
def resolve_path(x):
	return(str(pathlib2.Path(x).resolve(strict=False)))

##############
# containers #
##############

bbduk_container = "shub://TomHarrop/singularity-containers:bbmap_38.00"
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'
busco_container = "docker://ezlabgva/busco:v5.2.1_cv1"
salmon_container = "docker://combinelab/salmon:1.5.1"

#########
# rules #
#########

rule target:
	input:
		expand("output/FASTQC_untrimmed_short_reads/short_reads_{sample}_R{n}_fastqc.zip", sample=["01", "02"], n=[1,2]),
		expand("output/FASTQC_trimmed_short_reads/short_reads_{sample}_R{n}_fastqc.zip", sample=["01", "02"], n=[1,2]),
		"output/bbstats_hybrid/bbstats.out",
		"output/busco/hybrid_assembly/short_summary.specific.insecta_odb10.hybrid_assembly.txt",
		"output/emapper/cricket_hybrid_emapper.emapper.annotations.zst",
		"output/salmon/hybrid_transcriptome_quant/quant.sf",
		"output/bbstats_long/bbstats_long.out",
		"BUSCO_diagram/busco_combined_figure.png"


#cat short read files togetehr
rule cat_short_reads:
	output:
		s01_R1 = "short_reads/short_reads_01_R1.fastq.gz", 
		s01_R2 = "short_reads/short_reads_01_R2.fastq.gz",
		s02_R1 = "short_reads/short_reads_02_R1.fastq.gz",
		s02_R2 = "short_reads/short_reads_02_R2.fastq.gz"
	shell:
		"cat short_reads/*/*-6267-43-49-01_S24_L00*_R1_001.fastq.gz > {output.s01_R1} & "
		"cat short_reads/*/*-6267-43-49-01_S24_L00*_R2_001.fastq.gz > {output.s01_R2} & "
		"cat short_reads/*/*-6267-43-49-02_S25_L00*_R1_001.fastq.gz > {output.s02_R1} & "
		"cat short_reads/*/*-6267-43-49-02_S25_L00*_R2_001.fastq.gz > {output.s01_R2} & "
		"wait"

#FASTQC quality check 
rule FASTQC_untrimmed_short_reads:
	input:
		expand('short_reads/short_reads_{sample}_R{n}.fastq.gz', sample=["01", "02"], n=[1,2])
	output:
		expand("output/FASTQC_untrimmed_short_reads/short_reads_{sample}_R{n}_fastqc.zip", sample=["01", "02"], n=[1,2])
	params:
		outdir = directory("output/FASTQC_untrimmed_short_reads")
	singularity:
		fastqc_container
	shell:
		"mkdir -p {params.outdir} ; "
		"fastqc --outdir {params.outdir} {input}"

rule bbduk_trim:
	input:
		R1 = "short_reads/short_reads_{sample}_R1.fastq.gz",
		R2 = "short_reads/short_reads_{sample}_R2.fastq.gz"
	output:
		R1 = "output/bbduk_trim/short_reads_{sample}_R1.fq.gz",
		R2 = "output/bbduk_trim/short_reads_{sample}_R2.fq.gz"
	params:
		adapters = "/adapters.fa"
	log:
		"output/log/bbduk_{sample}.log"
	singularity:
		bbduk_container
	shell:
		"bbduk.sh in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 &> {log}"

rule FASTQC_trimmed_short_reads:
	input:
		expand('output/bbduk_trim/short_reads_{sample}_R{n}.fq.gz', sample=["01", "02"], n=[1,2])
	output:
		expand("output/FASTQC_trimmed_short_reads/short_reads_{sample}_R{n}_fastqc.zip", sample=["01", "02"], n=[1,2])
	params:
		outdir = directory("output/FASTQC_trimmed_short_reads")
	singularity:
		fastqc_container
	shell:
		"mkdir -p {params.outdir} ; "
		"fastqc --outdir {params.outdir} {input}"

#########################################
## HYBRID TRANSCRIPTOME WITH RNASPADES ##
#########################################

#running rnaSPADES
rule run_hybrid_transcriptome:
	input:
		s1 = "output/bbduk_trim/short_reads_01_R1.fq.gz",
		s2 = "output/bbduk_trim/short_reads_01_R2.fq.gz",
		nanopore = "data/reads.scrubb.fastq.gz"
	output:
		"output/done_rnaspades/transcripts.fasta"
	params:
		wd = "output/done_rnaspades"
	threads:
		30
	log:
		"output/log/rnaspades.log"
	shell:
		"bin/SPAdes-3.15.3-Linux/bin/rnaspades.py -1 {input.s1} -2 {input.s2} --nanopore {input.nanopore} --ss rf -t {threads} -o {params.wd} &> {log}"

rule bbstats:
	input:
		"output/done_rnaspades/transcripts.fasta"
	output:
		"output/bbstats_hybrid/bbstats.out"
	log:
		"output/log/bbstats.log"
	singularity:
		bbduk_container
	shell:
		"stats.sh "
		"in={input} "
		"out={output} "
		"&> {log}"

rule busco_hybrid:
	input:
		transcripts = "output/done_rnaspades/transcripts.fasta"
	output:
		"output/busco/hybrid_assembly/short_summary.specific.insecta_odb10.hybrid_assembly.txt"
	log:
		str(pathlib2.Path(resolve_path('output/logs/'),
							'busco_hybrid.log'))
	params:
		wd = "output/busco",
		outdir = "hybrid_assembly",
		hybrid_transcripts = lambda wildcards, input: resolve_path(input.transcripts)
	singularity:
		busco_container
	threads:
		25
	shell:
		"cd {params.wd} || exit 1 ; "
		"busco "
		"--in {params.hybrid_transcripts} "
		"--out {params.outdir} "
		"--lineage insecta_odb10 "
		"--cpu {threads} "
		"--mode transcriptome "
		"-f "
		"&> {log}"

rule cricket_hybrid_emapper:
	input:
		hybrid_transcriptome = "output/done_rnaspades/transcripts.fasta"
	output:
		"output/emapper/cricket_hybrid_emapper.emapper.annotations.zst"
	params:
		wd = "output/emapper",
		prefix = "cricket_hybrid_emapper"
	threads:
		64
	log:
		"output/log/emapper.log"
	shell:
		"python3 /Volumes/archive/deardenlab/guhlin/software/eggnog-mapper-2.0.1/emapper.py "
		"--cpu {threads} "
		"-i {input.hybrid_transcriptome} "
		"--itype CDS "
		"-o {params.prefix} "
		"--output_dir {params.wd} "
		"&> {log}"

rule salmon_index:
	input:
		hybrid_transcriptome = "output/done_rnaspades/transcripts.fasta"
	output:
		"output/salmon/transcripts_index/refseq.bin"
	params:
		outdir = "output/salmon/transcripts_index"
	threads:
		20
	singularity:
		salmon_container
	log:
		"output/logs/salmon_index.log"
	shell:
		"salmon index "
		"-t {input.hybrid_transcriptome} "
		"-i {params.outdir} "
		"-p {threads} "
		"&> {log}"
		
rule salmon_back:
	input:
		index_output = "output/salmon/transcripts_index/refseq.bin",
		trimmed_r1 = "output/bbduk_trim/short_reads_01_R1.fq.gz",
		trimmed_r2 = "output/bbduk_trim/short_reads_01_R2.fq.gz"
	output:
		quant = "output/salmon/hybrid_transcriptome_quant/quant.sf"
	params:
		index_outdir = "output/salmon/transcripts_index",
		outdir = "output/salmon/hybrid_transcriptome_quant"
	threads:
		20
	singularity:
		salmon_container
	log:
		"output/salmon/salmon_logs/salmon_quant.log"
	shell:
		"salmon quant "
		"-i {params.index_outdir} "
		"-l ISR "
		"-1 {input.trimmed_r1} "
		"-2 {input.trimmed_r2} "
		"-o {params.outdir} "
		"-p {threads} "
		"&> {log}"


######################################
## LONG READ TRANSCRIPTOME ANALYSIS ##
######################################

rule bbstats_long:
	input:
		"data/rnabloom.transcripts.fa"
	output:
		"output/bbstats_long/bbstats_long.out"
	log:
		"output/log/bbstats_long.log"
	singularity:
		bbduk_container
	shell:
		"stats.sh "
		"in={input} "
		"out={output} "
		"&> {log}"


#make plot of busco summaries
rule plot_busco:
	input:
		ss = expand("BUSCO_diagram/short_summary.specific.insecta_odb10.{type}.txt",
			type=["hybrid_assembly", "long_assembly"])
	output:
		busco_plot = "BUSCO_diagram/busco_combined_figure.png"
	params:
		ss_dir = "BUSCO_diagram/"
	singularity:
		busco_container
	threads:
		25
	shell:
		"python3 scripts/Busco_plot.py "
		'-wd {params.ss_dir}'

