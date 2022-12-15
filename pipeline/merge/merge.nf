#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.gasm_vcf = '/home/egarcia/data/igv_workfiles/vcf/DBA_2J-gasm.vcf'
params.pileup_vcf = '/home/egarcia/data/igv_workfiles/vcf/DBA_2J.snps-indels.pileup-cram.vcf'
params.results = "./results"
params.merge_script = "./bin/merge.py"

process MERGE {

	maxForks 5
	publishDir file(params.results + '/merged'), mode: "copy"

    input:
        tuple val(strain), file(gasm_vcf), file(pileup_vcf)

    output:
		tuple val(strain), file("*.vcf")

	script:
		merge = file(params.merge_script)

    shell:
    '''
	bcftools view -h !{gasm_vcf} > "!{strain}-merged.vcf"

	sed -i '$ d' "!{strain}-merged.vcf"

	echo "##INFO=<ID=CTYPE,Number=1,Type=String,Description="Calculated SVTYPE">" >> "!{strain}-merged.vcf"
	echo "##INFO=<ID=CLEN,Number=1,Type=Integer,Description="Calculated SVLEN">" >> "!{strain}-merged.vcf"
	echo "##INFO=<ID=CEND,Number=1,Type=Integer,Description="Calculated END">" >> "!{strain}-merged.vcf"
	echo "##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">" >> "!{strain}-merged.vcf"
	echo "##INFO=<ID=GASM,Type=String,Description="Coordinates for GASM.">" >> "!{strain}-merged.vcf"
	echo "##INFO=<ID=PILEUP,Type=String,Description="Coordinates for PILEUP.">" >> "!{strain}-merged.vcf"

	echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample" >> "!{strain}-merged.vcf"

	mkdir "!{strain}_tmp"
	export TMPDIR="!{strain}_tmp"

    python !{merge} --gasm_file=!{gasm_vcf} \
        --pileup_file=!{pileup_vcf} --out=!{strain}-merged.vcf.tmp

    cat !{strain}-merged.vcf.tmp >> "!{strain}-merged.vcf"
    '''

}

process INDELS {

	maxForks 5
	publishDir file(params.results + '/merged'), mode: "copy"

	input:
		tuple val(strain), file(merged_vcf)

	output:
		file "*.bed"
		file "*.csv"

	shell:
	'''
	bcftools query -i'CTYPE="DEL"' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} > !{strain}-merged.DEL.bed
	bcftools query -i'CTYPE="INS"' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} > !{strain}-merged.INS.bed

	bcftools query -i'CTYPE="INS" | CTYPE="DEL"' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} > !{strain}-merged.bed

	DEL_0_10="$(bcftools query -i'CTYPE="DEL" & abs(CLEN)<10' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	DEL_10_30="$(bcftools query -i'CTYPE="DEL" & abs(CLEN)>=10 & abs(CLEN)<30' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	DEL_30_50="$(bcftools query -i'CTYPE="DEL" & abs(CLEN)>=30 & abs(CLEN)<50' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	DEL_50_100="$(bcftools query -i'CTYPE="DEL" & abs(CLEN)>=50 & abs(CLEN)<100' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	DEL_100_1000="$(bcftools query -i'CTYPE="DEL" & abs(CLEN)>=100 & abs(CLEN)<1000' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	DEL_1000_10000="$(bcftools query -i'CTYPE="DEL" & abs(CLEN)>=1000 & abs(CLEN)<10000' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	DEL_10000_100000="$(bcftools query -i'CTYPE="DEL" & abs(CLEN)>=10000 & abs(CLEN)<100000' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	DEL_100000_MORE="$(bcftools query -i'CTYPE="DEL" & abs(CLEN)>=100000' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"

	INS_0_10="$(bcftools query -i'CTYPE="INS" & abs(CLEN)<10' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	INS_10_30="$(bcftools query -i'CTYPE="INS" & abs(CLEN)>=10 & abs(CLEN)<30' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	INS_30_50="$(bcftools query -i'CTYPE="INS" & abs(CLEN)>=30 & abs(CLEN)<50' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	INS_50_100="$(bcftools query -i'CTYPE="INS" & abs(CLEN)>=50 & abs(CLEN)<100' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	INS_100_1000="$(bcftools query -i'CTYPE="INS" & abs(CLEN)>=100 & abs(CLEN)<1000' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	INS_1000_10000="$(bcftools query -i'CTYPE="INS" & abs(CLEN)>=1000 & abs(CLEN)<10000' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	INS_10000_100000="$(bcftools query -i'CTYPE="INS" & abs(CLEN)>=10000 & abs(CLEN)<100000' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"
	INS_100000_MORE="$(bcftools query -i'CTYPE="INS" & abs(CLEN)>100000' -f'%CHROM\t%POS\t%CEND\t%CLEN\n' !{merged_vcf} | wc -l)"

	echo -e "!{strain},DEL,$DEL_0_10,$DEL_10_30,$DEL_30_50,$DEL_50_100,$DEL_100_1000,$DEL_1000_10000,$DEL_10000_100000,$DEL_100000_MORE" >> !{strain}-stats.csv
	echo -e "!{strain},INS,$INS_0_10,$INS_10_30,$INS_30_50,$INS_50_100,$INS_100_1000,$INS_1000_10000,$INS_10000_100000,$INS_100000_MORE" >> !{strain}-stats.csv
	'''

}

workflow {

	Channel.fromPath(params.gasm_vcf).map { file ->
		def strain = file.name.tokenize('-').get(0)
		return tuple(strain, file)
	}.set { gasm_files }

	Channel.fromPath(params.pileup_vcf).map { file ->
		def strain = file.name.tokenize('.').get(0)
		return tuple(strain, file)
	}.set { pileup_files }

	gasm_files.join(pileup_files, by: 0).set {merge_input}


	merged = MERGE( merge_input )
	INDELS( merged )
}

