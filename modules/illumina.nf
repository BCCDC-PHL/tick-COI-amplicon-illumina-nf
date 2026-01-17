process fastp {

    tag { sampleName }
    
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_fastp.json", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("${sampleName}_fastp.{json,html}"), emit: json
    tuple val(sampleName), path("${sampleName}_trimmed_R1.fastq.gz"), path("${sampleName}_trimmed_R2.fastq.gz"), emit: fastp_trimmed_reads
    tuple val(sampleName), path("${sampleName}_fastp_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastp\\n"  >> ${sampleName}_fastp_provenance.yml
    printf -- "  tools:\\n"               >> ${sampleName}_fastp_provenance.yml
    printf -- "    - tool_name: fastp\\n" >> ${sampleName}_fastp_provenance.yml
    printf -- "      tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sampleName}_fastp_provenance.yml
    
    fastp \
    --cut_tail \
    --trim_poly_g \
	-i ${forward} \
	-I ${reverse} \
    --detect_adapter_for_pe \
	-o ${sampleName}_trimmed_R1.fastq.gz \
	-O ${sampleName}_trimmed_R2.fastq.gz \
    --json ${sampleName}_fastp.json \
    --html ${sampleName}_fastp.html
        
    """
}

process performHostFilter {

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_hostfiltered_R*.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_performHostFilter_provenance.yml", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("${sampleName}_hostfiltered_R1.fastq.gz"), path("${sampleName}_hostfiltered_R2.fastq.gz"), emit: fastqPairs
    tuple val(sampleName), path("${sampleName}_performHostFilter_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: performHostFilter\\n"                                           >> ${sampleName}_performHostFilter_provenance.yml
    printf -- "  tools:\\n"                                                                    >> ${sampleName}_performHostFilter_provenance.yml
    printf -- "    - tool_name: bwa mem\\n"                                                    >> ${sampleName}_performHostFilter_provenance.yml
    printf -- "      tool_version: \$(bwa 2>&1 | grep "Version: " | cut -d ' ' -f 2)\\n"       >> ${sampleName}_performHostFilter_provenance.yml
    printf -- "      parameters:\\n"                                                           >> ${sampleName}_performHostFilter_provenance.yml
    printf -- "        - parameter: -t\\n"                                                   >> ${sampleName}_performHostFilter_provenance.yml
    printf -- "          value: ${task.cpus}\\n"                                             >> ${sampleName}_performHostFilter_provenance.yml

    bwa mem -t ${task.cpus} ${params.composite_ref} ${forward} ${reverse} | \
      filter_non_human_reads.py -c ${params.viral_contig_name} > ${sampleName}.viral_and_nonmapping_reads.bam
    samtools sort -@ ${task.cpus} -n ${sampleName}.viral_and_nonmapping_reads.bam | \
      samtools fastq -1 ${sampleName}_hostfiltered_R1.fastq.gz -2 ${sampleName}_hostfiltered_R2.fastq.gz -s ${sampleName}_singletons.fastq.gz -
    """
}

process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_readTrimming*", mode: 'copy'

    cpus 1

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), optional: true, emit: trimmedReads
    tuple val(sampleName), path("${sampleName}_readTrimming_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: readTrimming\\n"                                                                     >> ${sampleName}_readTrimming_provenance.yml
    printf -- "  tools:\\n"                                                                                         >> ${sampleName}_readTrimming_provenance.yml
    printf -- "    - tool_name: trim_galore\\n"                                                                     >> ${sampleName}_readTrimming_provenance.yml
    printf -- "      tool_version: \$(trim_galore --version | sed -n '4p' | sed 's/version //' | tr -d '[:space:]')\\n"                 >> ${sampleName}_readTrimming_provenance.yml

    trim_galore --paired $forward $reverse
    """
}

process filterResidualAdapters {
    /**
    * Discard reads that contain residual adapter sequences that indicate trimming may have failed
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output untrim_filter_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*{1,2}_posttrim_filter.fq.gz', mode: 'copy'

    cpus 1

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*1_posttrim_filter.fq.gz"), path("*2_posttrim_filter.fq.gz")

    script:
    """
    filter_residual_adapters.py --input_R1 $forward --input_R2 $reverse
    """
}

process kraken2Reports {
    /**Run Kraken2
    */
    tag { sampleName }
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: "${sampleName}_kraken2_report.txt", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: "${sampleName}_kraken2Reports*", mode: 'copy'


    input:
    tuple val(sampleName), path(forward), path(reverse), path(kraken2_db)

    output:
    path("${sampleName}_kraken2_report.txt")
    tuple val(sampleName), path("${sampleName}_kraken2_output.tsv"), path("${sampleName}_kraken2_report.txt"), optional:true, emit: report
    tuple val(sampleName), path("${sampleName}_kraken2Reports_provenance.yml"), emit: provenance


    script:
    """
    printf -- "- process_name: kraken2Reports\\n"                                              >> ${sampleName}_kraken2Reports_provenance.yml
    printf -- "  tools:\\n"                                                                    >> ${sampleName}_kraken2Reports_provenance.yml
    printf -- "    - tool_name: kraken2\\n"                                                    >> ${sampleName}_kraken2Reports_provenance.yml
    printf -- "      tool_version: \$(kraken2 --version | cut -d ' ' -f 3 | head -n 1)\\n"     >> ${sampleName}_kraken2Reports_provenance.yml
    printf -- "      parameters:\\n"                                                           >> ${sampleName}_kraken2Reports_provenance.yml
    printf -- "        - parameter: --db\\n"                                               >> ${sampleName}_kraken2Reports_provenance.yml
    printf -- "          value: \$(readlink ${kraken2_db})\\n"                                 >> ${sampleName}_kraken2Reports_provenance.yml
    printf -- "        - parameter: --threads\\n"                                                   >> ${sampleName}_kraken2Reports_provenance.yml
    printf -- "          value: ${task.cpus}\\n"                                               >> ${sampleName}_kraken2Reports_provenance.yml
    
    kraken2 \
      --db ${kraken2_db} \
      --threads ${task.cpus} \
      --output ${sampleName}_kraken2_output.tsv \
      --report ${sampleName}_kraken2_report.txt \
      --paired ${forward} ${reverse}
    """
}

process removeBacterialReads {
    /** extract all reads NOT classified as Bacteria (taxid 2) or any classification in the Bacteria subtree
    */
  
    tag { sampleName }
    errorStrategy 'ignore'

    publishDir "${params.outdir}/coiIllumina_sequenceAnalysis_kraken2Reports/${sampleName}", pattern: "kraken2_report_summary.csv", mode: 'copy'
    publishDir "${params.outdir}/coiIllumina_sequenceAnalysis_kraken2Reports/${sampleName}", pattern: "${sampleName}_removeBacterialReads_provenance.yml", mode: 'copy'

    cpus 1

    input:
    tuple val(sampleName), path(forward), path(reverse), path(kraken2_output), path(kraken2_report)

    output:
    tuple val(sampleName), path("${sampleName}_hostfiltered_R1_val_1_posttrim_filter_unclassified.fq.gz"), path("${sampleName}_hostfiltered_R2_val_2_posttrim_filter_unclassified.fq.gz"), emit: filtered
    path("kraken2_report_summary.csv"), emit: kraken2_report
    tuple val(sampleName), path("${sampleName}_removeBacterialReads_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: removeBacterialReads\\n"                                        >> ${sampleName}_removeBacterialReads_provenance.yml
    printf -- "  tools:\\n"                                                                    >> ${sampleName}_removeBacterialReads_provenance.yml
    printf -- "    - tool_name: extract_kraken_reads.py\\n"                                    >> ${sampleName}_removeBacterialReads_provenance.yml
    printf -- "      parameters:\\n"                                                           >> ${sampleName}_removeBacterialReads_provenance.yml
    printf -- "        - parameter: --taxid\\n"                                                  >> ${sampleName}_removeBacterialReads_provenance.yml
    printf -- "          value: 2\\n"                                                          >> ${sampleName}_removeBacterialReads_provenance.yml
    printf -- "        - parameter: --exclude\\n"                                                >> ${sampleName}_removeBacterialReads_provenance.yml
    printf -- "        - parameter: --include-children\\n"                                       >> ${sampleName}_removeBacterialReads_provenance.yml
    printf -- "        - parameter: --fastq-output\\n"                                           >> ${sampleName}_removeBacterialReads_provenance.yml

    extract_kraken_reads.py \
      -k ${kraken2_output} \
      -r ${kraken2_report} \
      -1 ${forward} \
      -2 ${reverse} \
      --taxid 2 \
      --exclude \
      --include-children \
      --fastq-output \
      -o ${sampleName}_hostfiltered_R1_val_1_posttrim_filter_unclassified.fq.gz \
      -o2 ${sampleName}_hostfiltered_R2_val_2_posttrim_filter_unclassified.fq.gz

    kraken2_report.py \
      --report ${kraken2_report} \
      --output_file kraken2_report_summary.csv \
      --samplename ${sampleName} \
      --prefix ${params.prefix}
    """
}

process subsample {
        /**
    * subsample to a x number of reads (default is 1,000)
    * 
    */
    tag { sampleName }
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*subsampled.fastq.gz', mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_subsample_provenance.yml", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("${sampleName}.R1.subsampled.fastq.gz"), path("${sampleName}.R2.subsampled.fastq.gz"), emit: subsample_out
    tuple val(sampleName), path("${sampleName}_subsample_provenance.yml"), emit: provenance


    script:
    """
    printf -- "- process_name: subsample\\n"                                                   >> ${sampleName}_subsample_provenance.yml
    printf -- "  tools:\\n"                                                                    >> ${sampleName}_subsample_provenance.yml
    printf -- "    - tool_name: rasusa\\n"                                                     >> ${sampleName}_subsample_provenance.yml
    printf -- "      tool_version: \$(rasusa --version | sed 's/rasusa //')\\n"                >> ${sampleName}_subsample_provenance.yml
    printf -- "      parameters:\\n"                                                           >> ${sampleName}_subsample_provenance.yml
    printf -- "        - parameter: --num\\n"                                                    >> ${sampleName}_subsample_provenance.yml
    printf -- "          value: ${params.subsampleReads}\\n"                                   >> ${sampleName}_subsample_provenance.yml
    printf -- "        - parameter: --seed\\n"                                                   >> ${sampleName}_subsample_provenance.yml
    printf -- "          value: 42\\n"                                                         >> ${sampleName}_subsample_provenance.yml

    rasusa --seed 42 --num ${params.subsampleReads} -i ${forward} -i ${reverse} -o ${sampleName}.R1.subsampled.fastq.gz -o ${sampleName}.R2.subsampled.fastq.gz
    """

}

process denovoAssembly {
        /**
    * assembly using SPAdes  
    * 
    */

    tag { sampleName }
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: 'results*', mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: '*.fasta', mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: '*.json', mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: "${sampleName}_denovoAssembly_provenance.yml", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    path("results*")
    path("*.fasta")
    path("*.json")
    tuple val(sampleName), path("${sampleName}_*_filtered.fasta"), optional: true,  emit: scaffolds
    tuple val(sampleName), path("${sampleName}_denovoAssembly_provenance.yml"), emit: provenance


    script:
    """
    printf -- "- process_name: denovoAssembly\\n"                                              >> ${sampleName}_denovoAssembly_provenance.yml
    printf -- "  tools:\\n"                                                                    >> ${sampleName}_denovoAssembly_provenance.yml
    printf -- "    - tool_name: spades.py\\n"                                                  >> ${sampleName}_denovoAssembly_provenance.yml
    printf -- "      tool_version: \$(spades.py --version | cut -d ' ' -f 4)\\n"               >> ${sampleName}_denovoAssembly_provenance.yml
    printf -- "      parameters:\\n"                                                           >> ${sampleName}_denovoAssembly_provenance.yml
    printf -- "        - parameter: -t\\n"                                                   >> ${sampleName}_denovoAssembly_provenance.yml
    printf -- "          value: ${task.cpus}\\n"                                               >> ${sampleName}_denovoAssembly_provenance.yml
    printf -- "        - parameter: --isolate\\n"                                                >> ${sampleName}_denovoAssembly_provenance.yml

    spades.py -1 ${forward} -2 ${reverse} -o results/${sampleName}/ -t ${task.cpus}

        if [ -e results/${sampleName}/scaffolds.fasta ]; then
            cp results/${sampleName}/scaffolds.fasta ${sampleName}.scaffolds.fasta
            filter_seqs.py --input_file ${sampleName}.scaffolds.fasta --output_file ${sampleName}.fasta --length ${params.length} --output_file_rename ${sampleName}_filtered.fasta --samplename ${sampleName}
            
            if [ -e ${sampleName}_NODE_*_filtered.fasta ]; then
                touch complete_assembly.json
            else
                rm *_filtered.fasta
                spades.py -1 ${forward} -2 ${reverse} -o results_retry/${sampleName}/ -t ${task.cpus} --isolate

                if [ -e results_retry/${sampleName}/scaffolds.fasta ]; then
                    cp results_retry/${sampleName}/scaffolds.fasta ${sampleName}.scaffolds.fasta
                    filter_seqs.py --input_file ${sampleName}.scaffolds.fasta --output_file ${sampleName}.fasta --length ${params.length} --output_file_rename ${sampleName}_filtered.fasta --samplename ${sampleName}
                    if [ -e ${sampleName}_NODE_*_filtered.fasta ]; then
                        touch complete_assembly.json
                    else
                        touch incomplete_assembly.json
                    fi
                fi
            fi
        else
            touch empty.json
        fi
    """

}

process alignConsensusToReference {
    /**
    * Aligns consensus sequence against reference using mafft. Uses the --keeplength
    * flag to guarantee that all alignments remain the same length as the reference.
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: "${sampleName}.with_ref.alignment.fa", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/consensus", pattern: "${sampleName}.consensus.fasta", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: "${sampleName}_alignConsensusToReference_provenance.yml", mode: 'copy'

    input:
    tuple val(sampleName), path(consensus), path(reference)

    output:
    tuple val(sampleName), path("${sampleName}.with_ref.alignment.fa")
    tuple val(sampleName), path("${sampleName}.consensus.fasta"), emit: consensus
    tuple val(sampleName), path("${sampleName}_alignConsensusToReference_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: alignConsensusToReference\\n"                                   >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "  tools:\\n"                                                                    >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "    - tool_name: mafft\\n"                                                      >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "      tool_version: \$(mafft --version 2>&1)\\n"                                     >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "      parameters:\\n"                                                           >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "        - parameter: --adjustdirectionaccurately\\n"                              >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "        - parameter: --preservecase\\n"                                           >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "        - parameter: --keeplength\\n"                                             >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "        - parameter: --add\\n"                                                    >> ${sampleName}_alignConsensusToReference_provenance.yml
    printf -- "          value: \$(readlink ${reference})\\n"                                  >> ${sampleName}_alignConsensusToReference_provenance.yml

    mafft \
      --adjustdirectionaccurately \
      --preservecase \
      --keeplength \
      --add \
      ${consensus} \
      ${reference} \
      > ${sampleName}.with_ref.alignment.fa
    extract_seqs.py --alignment_file ${sampleName}.with_ref.alignment.fa --samplename ${sampleName}
    mv *${sampleName}*.consensus.fasta ${sampleName}.consensus.fasta
    """

}

process blastSpeciesID {
    /** Perform BLASTN to determine Ixodes/Dermacentor tick species ID **/

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_blastn.csv", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_blastSpeciesID_provenance.yml", mode: 'copy'
    
    input:
    tuple val(sampleName), path(consensus), path(blastndb)

    output:
    tuple val(sampleName), path("${sampleName}_blastn.csv"), emit: blastn
    tuple val(sampleName), path("${sampleName}_blastSpeciesID_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: blastSpeciesID\\n"                                                 >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "  tools:\\n"                                                                       >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "    - tool_name: blastn\\n"                                                        >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "      tool_version: \$(blastn -version | cut -d ' ' -f 2 | head -n 1)\\n"          >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "      parameters:\\n"                                                              >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "        - parameter: -db\\n"                                                        >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "          value: \$(readlink ${blastndb})\\n"                                      >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "          sha256: \$(shasum -a 256 ${blastndb}/*.fasta | cut -d ' ' -f 1)\\n"              >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "        - parameter: -outfmt\\n"                                                    >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "          value: 10 qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore\\n"     >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "        - parameter: -max_target_seqs\\n"                                           >> ${sampleName}_blastSpeciesID_provenance.yml
    printf -- "          value: ${params.max_target_seqs_blast}\\n"                               >> ${sampleName}_blastSpeciesID_provenance.yml

    grep "^>" ${consensus} > ${sampleName}_consensus_temp.fa && \
    grep -v "^>" ${consensus} | tr -d "-" | tr -d "\n" | fold -w 70 >> ${sampleName}_consensus_temp.fa

    if [ \$(echo \$(grep -v '>' ${sampleName}_consensus_temp.fa | tr -d '\n' | wc -c)) -gt 0 ]; then
    	blastn \
    	-db ${blastndb}/*.fasta \
    	-query ${sampleName}_consensus_temp.fa \
    	-outfmt "10 qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore" \
      	-max_target_seqs ${params.max_target_seqs_blast} \
    	-out ${sampleName}_blastn.csv
    	echo "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore" > header.csv
    	cat header.csv ${sampleName}_blastn.csv > temp.csv
    	mv temp.csv ${sampleName}_blastn.csv
    	rm header.csv ${sampleName}_consensus_temp.fa
    else
        touch ${sampleName}_blastn.csv
    fi
    """
}

process indexReferences {
    /**
    * Indexes the consensus sequences produced in alignConsensusToReference
    */

    tag { ref }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: "*", mode: 'copy'

    input:
    tuple val(sampleName), path(ref)

    output:
    tuple val(sampleName), path("${ref}.fa"), path("${ref}.fa.*"), emit: indexed_ref
    tuple val(sampleName), path("${sampleName}_indexReferences_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: indexReferences\\n"                                                 >> ${sampleName}_indexReferences_provenance.yml
    printf -- "  tools:\\n"                                                                        >> ${sampleName}_indexReferences_provenance.yml
    printf -- "    - tool_name: bwa\\n"                                                            >> ${sampleName}_indexReferences_provenance.yml
    printf -- "      tool_version: \$(bwa 2>&1 | grep "Version: " | cut -d ' ' -f 2)\\n"           >> ${sampleName}_indexReferences_provenance.yml
    printf -- "      parameters:\\n"                                                               >> ${sampleName}_indexReferences_provenance.yml
    printf -- "        - parameter: index\\n"                                                      >> ${sampleName}_indexReferences_provenance.yml

    ln -s ${ref} ${ref}.fa
    bwa index ${ref}.fa
    """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: "${sampleName}.sorted{.bam,.bam.bai}", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/${sampleName}", pattern: "${sampleName}_readMapping_provenance.yml", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse), path(ref), path("*")

    output:
    tuple val(sampleName), path("${sampleName}.sorted.bam"), path("${sampleName}.sorted.bam.bai"), emit: bamfiles
    tuple val(sampleName), path("${sampleName}_readMapping_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: readMapping\\n"                                                     >> ${sampleName}_readMapping_provenance.yml
    printf -- "  tools:\\n"                                                                        >> ${sampleName}_readMapping_provenance.yml
    printf -- "    - tool_name: bwa\\n"                                                            >> ${sampleName}_readMapping_provenance.yml
    printf -- "      tool_version: \$(bwa 2>&1 | grep "Version: " | cut -d ' ' -f 2)\\n"           >> ${sampleName}_readMapping_provenance.yml
    printf -- "      parameters:\\n"                                                               >> ${sampleName}_readMapping_provenance.yml
    printf -- "        - parameter: mem\\n"                                                        >> ${sampleName}_readMapping_provenance.yml
    printf -- "    - tool_name: samtools\\n"                                                       >> ${sampleName}_readMapping_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep "Version: " | cut -d ' ' -f 2)\\n"      >> ${sampleName}_readMapping_provenance.yml
    printf -- "      parameters:\\n"                                                               >> ${sampleName}_readMapping_provenance.yml
    printf -- "        - parameter: sort\\n"                                                       >> ${sampleName}_readMapping_provenance.yml
    printf -- "        - parameter: index\\n"                                                      >> ${sampleName}_readMapping_provenance.yml

    bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | \
    samtools sort -o ${sampleName}.sorted.bam
    samtools index ${sampleName}.sorted.bam
    """
}

process extractSampleIDs {
    /** extract SampleIDs that have failed to meet the 'PASS' check_metrics criteria, 
    'pass' >= 658 bp in size and >= 97% sequence identity
    'review' >= 658 bp in size and <97% sequence identity
    'fail' <658 bp in size
     **/
    tag { params.prefix }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.prefix}_review_failed_samples.txt", mode: 'copy'
    
    input:
    path(summary_file)

    output:
    path("${params.prefix}_review_failed_samples.txt"), emit: extractSamples

    script:
    """
    extractSampleIDs.py --summary_file ${summary_file} --summary_outfile ${params.prefix}_review_failed_samples.txt
    """
}

process ncbiBlast {
    /** Perform BLASTN using the NCBI database to determine other non-Ixodes/Dermacentor species ID **/

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_blastn.tsv", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_ncbiBlast_provenance.yml", mode: 'copy'

    input:
    tuple val(sampleName), path(consensus)

    output:
    tuple val(sampleName), path("${sampleName}_blastn.tsv"), emit: blast_ncbi
    tuple val(sampleName), path("${sampleName}_ncbiBlast_provenance.yml"), emit: provenance


    script:
    """
    printf -- "- process_name: ncbiBlast\\n"                                                      >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "  tools:\\n"                                                                       >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "    - tool_name: blastn\\n"                                                        >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "      tool_version: \$(blastn -version | cut -d ' ' -f 2 | head -n 1)\\n"          >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "      parameters:\\n"                                                              >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "        - parameter: -db\\n"                                                        >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "          value: \$(readlink ${params.ncbi_db})\\n"                                >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "        - parameter: -outfmt\\n"                                                    >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "          value: 6 qseqid staxids sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue bitscore\\n"     >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "        - parameter: -max_target_seqs\\n"                                           >> ${sampleName}_ncbiBlast_provenance.yml
    printf -- "          value: ${params.max_target_seqs_blast}\\n"                               >> ${sampleName}_ncbiBlast_provenance.yml

    grep "^>" ${consensus} > ${sampleName}_consensus_temp.fa && \
    grep -v "^>" ${consensus} | tr -d "-" | tr -d "\n" | fold -w 70 >> ${sampleName}_consensus_temp.fa

    export BLASTDB=${params.ncbi_db}
    if [ \$(echo \$(grep -v '>' ${sampleName}_consensus_temp.fa | tr -d '\n' | wc -c)) -gt 0 ]; then
        blastn \
        -db ${params.ncbi_db}/nt \
        -query ${sampleName}_consensus_temp.fa \
        -outfmt "6 qseqid staxids sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -max_target_seqs ${params.max_target_seqs_blast} \
        -out ${sampleName}_blastn.tsv
    else
        touch ${sampleName}_blastn.tsv
    fi
    """
}

