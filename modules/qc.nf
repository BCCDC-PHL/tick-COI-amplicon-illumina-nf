process makeQCCSV {
    tag { sampleName }

    cpus 4

    input:
    tuple val(sampleName), path(bam), path(bam_index), path(fasta), path(ref), path(primer_bed), path(primer_pairs), path(top_hits_blastn)

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv

    script:
    """
    qc.py --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta} --primer-bed ${primer_bed} --primer-pairs ${primer_pairs} --min-depth ${params.varMinDepth} --top_hits_blastn ${top_hits_blastn} --prefix ${params.prefix}
    """
}

process writeQCSummaryCSV {
    tag { params.prefix }

    executor 'local'

    input:
    val lines

    output:
    path("${params.prefix}.summary.COI.qc.csv")

    exec:
    summary_COI = file("${params.outdir}/${params.prefix}.summary.COI.qc.csv")
    summary_COI.withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }
    }

    allLines = summary_COI.readLines()
    header_row = allLines[0]
    data_rows = allLines[1..-1].sort { it.split(',')[1] }
    sortedLines = ([header_row] + data_rows).join('\n')
    summary_COI.write(sortedLines)
    summary_COI.copyTo("${task.workDir}/${params.prefix}.summary.COI.qc.csv")

}


process compileTopHits {
    tag { params.prefix }

    publishDir "${params.outdir}/", pattern: "${params.prefix}_blastn_data_dictionary.csv", mode: 'copy'

    input:
    tuple val(sampleName), path(blastn)

    output:
    path "${params.prefix}.${sampleName}.top_hits_blastn.csv", emit: topHits
    tuple val(sampleName), path("${params.prefix}.${sampleName}.top_hits_blastn.csv"), emit: topHitsCSV
    path "${params.prefix}_blastn_data_dictionary.csv"

    script: 
    """
    blastcsv.py --outfile ${params.prefix}.${sampleName}.top_hits_blastn.csv --dictionary ${params.prefix}_blastn_data_dictionary.csv --sample ${sampleName} --blastn ${blastn}
    """
}

process writeTopHitsCSV {
    tag { params.prefix }

    executor 'local'

    input:
    val lines

    output:
    path("${params.prefix}_top_hits_blastn.csv"), emit: top_hits_blastn

    exec:
    topHits = file("${params.outdir}/${params.prefix}_top_hits_blastn.csv")
    topHits.withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
    topHits.copyTo("${task.workDir}/${params.prefix}_top_hits_blastn.csv")


}


process reportNCBI {
    tag { params.prefix }

    input:
    tuple val(sampleName), path(blastn)

    output:
    path "top_hits_blastn_ncbi.csv", emit: ncbi_topHits

    script: 
    """
    blast_ncbi_report.py --report ${blastn} --output_file top_hits_blastn_ncbi.csv --sample ${sampleName} --prefix ${params.prefix}
    """
}

process failedSummary {
    /** Provides a list of all the failed samples and the reasons for failing them **/

    tag { sampleName }

    publishDir "${params.outdir}/", pattern: "${params.prefix}.summary.COI.failed.csv", mode: 'copy'

    input:
    path(kraken2report), path(summary_coi_csv)

    output:
    path("${params.prefix}.summary.COI.failed.csv")

    script:
    """
    summary_failed_samples.py --samplesheet ${params.samplesheet} --summary_report ${summary_coi_csv} --kraken2_report ${kraken2report} --prefix ${params.prefix}
    """
}
