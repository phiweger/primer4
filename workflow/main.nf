nextflow.enable.dsl = 2


include { design; pcr; pseudo; cat} from './processes'


workflow {
    ch = channel.fromPath(params.input)
                .splitCsv(header:true)
                .map{ row-> tuple(row.name, row.method, row.variant) }


    // We don't pass various files as channels bc/ they all have indices which
    // we'd have to otherwise pass explicitly.
    // https://github.com/nextflow-io/hack17-varcall/issues/1
    chh = ch.combine(channel.fromPath(params.settings))
            .combine([params.data])
            .view()


    design(chh)
    
    pcr(design.out.map { it -> [it[3]] }.collect().view())
    
    pseudo(
        design.out
              .map { it -> [it[0], it[1], it[2]] } 
              .combine(pcr.out).view())

    cat(pseudo.out.collect())

    // blast(
    //     design.out
    //           .map { it -> tuple(it[0], it[3]) }
    //           .combine([params.blastdb])
    // )
}
