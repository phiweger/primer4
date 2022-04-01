process design {
    publishDir "${params.results}/${name}", pattern: 'log.txt', mode: 'copy', overwrite: true
    publishDir "${params.results}/${name}", pattern: "${name}.json", mode: 'copy', overwrite: true
    echo true
    errorStrategy 'ignore'
    memory '8 GB'
    cpus 8

    input:
        tuple(
            val(name),
            val(method),
            val(variant),
            path(config),
            val(data))

    output:
        tuple(
            val(name),
            val(method),
            path("${name}.json"),
            path("${name}.tsv"),
            path('log.txt'))

    script:
    """
    design.py -n ${name} -p ${name} -q '${variant}' -m ${method} -c ${config} -d ${data} 2>&1 > log.txt

    # Not sure why nf won't save the log when using 2> log.txt
    # Note the ' around the variant; otherwise will fail bc/ eg
    # NM_032682.6:c.484C>T would write to file T
    """
}


process pcr {
    memory '8 GB'
    cpus 8
    container 'nanozoo/ispcr:33--2df9365'

    input:
        path(primers)
        // path(reference)

    output:
        path('ispcr.fna')

    shell:
    '''
    cat !{primers} > primers.tsv

    file=$(jq -r '.reference' !{params.settings})
    reference=!{params.data}/${file}
    
    isPcr $reference primers.tsv -maxSize=4000 -minPerfect=15 -minGood=15 -out=fa ispcr.fna
    '''
}


// process blast {
//     publishDir "${params.results}/${name}", mode: 'copy', overwrite: true
//     cpus 8

//     input:
//         tuple(val(name), path(primers), val(blastdb))

//     output:
//         tuple(val(name), path('blast.tsv'))

//     """
//     blastn -task blastn -num_threads ${task.cpus} -query ${primers} -db ${blastdb} -evalue 100 -dust no -word_size 7 -outfmt 6 > blast.tsv
//     """
// }


process pseudo {
    input:
        tuple(val(name), val(method), path(candidates), path(annealing))

    output:
        path("${name}.tsv")

    """
    pseudo.py -m ${method} -c ${candidates} -a ${annealing} -o ${name}.tsv
    """
}


process cat {
    publishDir "${params.results}", mode: 'copy', overwrite: true

    input:
        path(results)

    output:
        path("summary.tsv")

    """
    cat ${results} | sort -k5,5 -k2,2 > summary.tsv
    """
}
