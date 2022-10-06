#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reference = "./input/Mus_musculus.GRCm39.dna.toplevel.fa"
params.results = "./results"
params.bubble_beds = "./input/bubble_bed/*.bed"
params.chromosomes_dir = "./input/chromosomes"

process build_pangenome {

    publishDir file(params.results + '/support-files/minigraph'), mode: "copy"
    
    input:
        file reference
        file all_chromosomes

    output:
        file "*.gfa"

    shell:
    '''
    echo "Using !{task.cpus} CPUs: "
    minigraph -xggs -t${task.cpus} ${reference} ${all_chromosomes} > pangenome_graph.gfa
    '''
}

process get_bubbles { 

    publishDir file(params.results + '/support-files/minigraph'), mode: "copy"
    
    input:
        file chromosomes
        each genome_graph
    
    output:
        tuple val(strain), file("*.bubbles.bed")

    script:

    strain = chromosomes.name.tokenize('.').get(0)

    """
    minigraph -xasm -t!{task.cpus} --call !{genome_graph} !{chromosomes} > !{strain}.bubbles.bed
    """
}

process candidate_indels {

    publishDir file(params.results + '/support-files/candidate-regions'), mode: "copy"

    input:
        tuple val(strain), file(bubbles_bed)
    
    output:
        tuple val(strain), file("*.tsv")

    shell:
    '''
    cat !{bubbles_bed} |  awk -F"[\t:#]" 'BEGIN {OFS = "\t"; print "#ref_chr","ref_pos","ref_end","len1","|","strain_chr","strain_pos","strain_end","len2","|","short","long","a","b","|","ref_npos","ref_nend","ref_nlen","|","strain_npos","strain_nend","strain_nlen","|","info"} 
        {
            seg=5; href=1; hi=$3-$2; lo=$13-$12; 
            if(hi<lo){
                href=0; lo=$3-$2; hi=$13-$12;
            }; 

            a1=hi%seg; a=(hi-a1)/seg; 
            b1=((seg+2)*a-lo)%2; b=((seg+2)*a-lo-b1)/2; 

            if(href){
                ref_x=$2-a; ref_y=$3+a; 
                str_x=$12-b; str_y=$13+b;
            }else{
                ref_x=$2-b; ref_y=$3+b; 
                str_x=$12-a; str_y=$13+a;
            }; 
            if($6!="." && $3-$2!=$7 && $1 ~ /^[0-9]*$/) 
                print $1,$2,$3,$3-$2,"|",$11,$12,$13,$7,"|",lo,hi,a,b,"|",ref_x,ref_y,ref_y-ref_x,"|",str_x,str_y,str_y-str_x,"|",($13-$12)-($3-$2)":"$7-($3-$2)":"$6
        }' > candidates-!{strain}.tsv    
    '''   
}

process reference_regions {

    publishDir file(params.results + '/support-files/reference-regions'), mode: "copy"

    input:
        tuple val(strain), file(candidates)
        file reference

    output:
        tuple val(strain), file("*.fa")

    shell:
    '''
    cut -f1,16,17,18 !{candidates} > candidates.bed
    bedtools getfasta -fi !{reference} -bed candidates.bed -fo !{strain}.ref.segments.fa
    '''

}

process strain_regions {

    publishDir file(params.results + '/support-files/candidate-regions'), mode: "copy"

    input:
        tuple val(strain), file(candidates)
        file chromosomes_dir
        val mask_pattern

    output:
        tuple val(strain), file("*.fa")
    
    script:
        mask = mask_pattern.replace("{strain}", strain)

    shell:
    '''
    cut -f1,20,21,22 !{candidates} | awk  'BEGIN {OFS="\t"} {pos=$2; if(pos<0){pos=0}; if($1 ~ /^[0-9]*$/) print "!{mask}"$1,pos,$3,$4}' > candidates.bed
    bedtools getfasta -fi !{chromosomes_dir}/!{strain}.fasta -bed candidates.bed -fo !{strain}.str.segments.fa
    sed -i "s/!{mask}//g" !{strain}.str.segments.fa

    cut -f1,16,17,20,21 !{candidates} | awk  'BEGIN {OFS="\t"} {pos=$4; if(pos<0){pos=0}; if($1 ~ /^[0-9]*$/) print ">"$1":"$2"-"$3,">"$1":"pos"-"$5}' > index.txt

    for line in index.txt
    do
        pos_ref="$(echo $line | cut -f1)"
        pos_str="$(echo $line | cut -f2)"
        sed -i "s/"$pos_str"/"$pos_ref"/g" !{strain}.str.segments.fa
    done
    '''

}

process align_segments {

    memory '8 GB'
    maxForks 1
    publishDir file(params.results + '/support-files/aligned-regions'), mode: "copy"

    input:
        tuple val(strain), file(segments)
        file reference

    output:
        tuple val(strain), file("*.bam")

    shell:
    '''
    mkdir tmp
    minimap2 -a -R '@RG\\tID:!{strain}\\tSM:!{strain}' --MD -Y -t !{task.cpus} !{reference} !{segments} | samtools view -bS - | samtools sort -T ./tmp -o !{strain}.segments.sorted.bam -
    '''
}

process asm_call {

    publishDir file(params.results + '/support-files/variant-calls'), mode: "copy"

    input:
        tuple val(strain), file(aligned_regions)
        file reference

    output:
        tuple val(strain), file("*.vcf")

    shell:
    '''
    samtools index -@ !{task.cpus} !{aligned_regions}
    svim-asm haploid ./ !{aligned_regions} !{reference}
    mv variants.vcf !{strain}-gasm.vcf
    mv sv-lengths.png !{strain}-gasm.lengths.png
    '''

}

process bed_files {

    publishDir file(params.results + '/calls/'), mode: "copy"

    input: 
        tuple val(strain), file(vcf_file)
        each type
    
    output:
        tuple val(strain), val(type), file('*.bed')

    script:

    simple_name = vcf_file.name.replace(".vcf","")
    strain = simple_name.tokenize('-').get(0)

    """
    bcftools query -i"SVTYPE='!{type}'" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' !{vcf_file} | \
            awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4}' > "!{strain}-gasm.!{type}.bed"
    """
}


workflow {

    reference = file(params.reference)
    chromosomes_dir = file(params.chromosomes_dir)

    bubbles = Channel.fromPath(params.bubble_beds).map{ it->
        def strain = it.name.tokenize(".").get(0)
        return tuple(strain, it)
    }

    cadidates = candidate_indels(bubbles)
    regions_strain = strain_regions(cadidates, chromosomes_dir,'{strain}#1#chr')
    aligned = align_segments(regions_strain, reference)

    vcf_files = asm_call(aligned, reference)
    sv_types = Channel.from(['INS', 'DEL', 'INV', 'DUP'])

    bed_files(vcf_files, sv_types)
    
}