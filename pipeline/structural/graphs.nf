#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reference = "./input/Mus_musculus.GRCm39.dna.toplevel.chr1.fa"
params.reads = "./input/reads/long/*/*.gz"
params.results = "./results"
params.bubble_beds = "./input/bubble_bed/*.bed"

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

    publishDir file(params.results + '/support-files/graphs/'), mode: "copy"

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
            if($3-$2!=$7 && $1 ~ /^[0-9]*$/) 
                print $1,$2,$3,$3-$2,"|",$11,$12,$13,$7,"|",lo,hi,a,b,"|",ref_x,ref_y,ref_y-ref_x,"|",str_x,str_y,str_y-str_x,"|",($13-$12)-($3-$2)":"$7-($3-$2)":"$6
        }' > candidates.tsv    
    '''   
}

workflow {
    bubbles = Channel.fromPath(params.bubble_beds).map{ it->
        def strain = it.name.tokenize(".").get(0)
        return tuple(strain, it)
    }
    candidate_indels(bubbles).view()
}