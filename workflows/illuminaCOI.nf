#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { performHostFilter }         from '../modules/illumina.nf'
include { readTrimming }              from '../modules/illumina.nf'
include { filterResidualAdapters }    from '../modules/illumina.nf'
include { kraken2Reports  }           from '../modules/illumina.nf'
include { removeBacterialReads  }     from '../modules/illumina.nf'
include { subsample }                 from '../modules/illumina.nf'
include { denovoAssembly  }           from '../modules/illumina.nf'
include { alignConsensusToReference } from '../modules/illumina.nf'
include { blastSpeciesID }	      from '../modules/illumina.nf'
include { indexReferences }            from '../modules/illumina.nf'
include { readMapping }               from '../modules/illumina.nf'
include { makeQCCSV }         from '../modules/qc.nf'
include { writeQCSummaryCSV } from '../modules/qc.nf'
include { compileTopHits }    from '../modules/qc.nf'
include { writeTopHitsCSV }   from '../modules/qc.nf'
include { extractSampleIDs  } from '../modules/illumina.nf'
include { ncbiBlast  } from '../modules/illumina.nf'
include { reportNCBI } from '../modules/qc.nf'

workflow prepareReferenceFiles {

    // Get reference fasta (without primers)
    if (params.ref) {
      Channel.fromPath(params.ref)
              .set{ ch_refFasta }
    } else {
      ch_refFasta = Channel.empty() {
        error "Reference not found. Please specify a valid path for the ref parameter."
      }
    }

    // Get reference fasta (with primers)
    if (params.ref_with_primers) {
      Channel.fromPath(params.ref_with_primers)
              .set{ ch_refFasta_primers }
    } else {
      ch_refFasta_primers = Channel.empty() {
        error "Reference with primers not found. Please specify a valid path for the ref-with-primers parameter."
      }
    }

    // Get blastndb if db parameter is supplied
    if (params.db) {
	Channel.fromPath(new File(params.db).getAbsolutePath())
		.set{ ch_blastndb }
    } else {
	ch_blastndb = Channel.empty() {
		error "Blast database not found. Please specify a valid path for the db parameter."
	}
    }

    // Get the primers and bed file
  ch_primerPairs = Channel.fromPath(params.primer_pairs_tsv)

    if (params.bed) {
      ch_bedFile = Channel.fromPath(params.bed)
    } else {
      ch_bedFile = Channel.empty() {
        error "Bed file not found. Please specify a valid path for bed parameter."
      }
    }


    // Get kraken2 db if parameter is supplied
    if (params.kraken_db) {
	Channel.fromPath(new File(params.kraken_db).getAbsolutePath())
		.set{ ch_krakendb }
    } else {
	ch_krakendb = Channel.empty() {
		error "Kraken2 database not found. Please specify a valid path for the kraken_db parameter."
	}
    }


    emit:
      blastndb = ch_blastndb
      kraken_db = ch_krakendb
      primer_pairs = ch_primerPairs
      bedfile = ch_bedFile
      refFasta = ch_refFasta
      refFasta_primers = ch_refFasta_primers

}


workflow sequenceAnalysis {
    take:
      ch_filePairs
      ch_refFasta
      ch_refFasta_primers
      ch_bedFile
      ch_primerPairs
      ch_blastndb
      ch_krakendb


    main:

      performHostFilter(ch_filePairs)

      readTrimming(performHostFilter.out.fastqPairs)

      filterResidualAdapters(readTrimming.out.trimmedReads)

      kraken2Reports(filterResidualAdapters.out.combine(ch_krakendb))

      removeBacterialReads(filterResidualAdapters.out.join(kraken2Reports.out.report))

      subsample(removeBacterialReads.out.filtered)

      denovoAssembly(subsample.out)

      alignConsensusToReference(denovoAssembly.out.scaffolds.combine(ch_refFasta))

      blastSpeciesID(alignConsensusToReference.out.consensus.combine(ch_blastndb))

      compileTopHits(blastSpeciesID.out.blastn)

      compileTopHits.out.topHits.splitCsv()
                       .unique()
                       .set { blastnTopHits }

      writeTopHitsCSV(blastnTopHits.toList())

      indexReferences(alignConsensusToReference.out.consensus)

      readMapping(subsample.out.join(indexReferences.out.indexed_ref))

      makeQCCSV(readMapping.out.bamfiles.join(alignConsensusToReference.out.consensus, by: 0)
                                   .combine(ch_refFasta_primers)
				   .combine(ch_bedFile)
				   .combine(ch_primerPairs)
				   .join(compileTopHits.out.topHitsCSV))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .set { qc }

      writeQCSummaryCSV(qc.toList())

      removeBacterialReads.out.kraken2_report.collectFile(keepHeader: true, sort: { it.text }, name: "${params.prefix}.kraken2_report_summary.csv", storeDir: "${params.outdir}")
      extractSampleIDs(writeQCSummaryCSV.out)

      extractSampleIDs.out.extractSamples.splitText()
          .map { it.trim() }
          .map { [it,] }
          .set { extract_samples_ch }

      ncbiBlast(alignConsensusToReference.out.consensus.join(extract_samples_ch))

      reportNCBI(ncbiBlast.out.blast_ncbi)
      
      reportNCBI.out.ncbi_topHits.collectFile(keepHeader: true, sort: { it.text }, name: "${params.prefix}_top_hits_blastn_ncbi.csv", storeDir: "${params.outdir}")

}

workflow coiIllumina {
    take:
      ch_filePairs
      ch_blastndb
      ch_krakendb


    main:
      prepareReferenceFiles()
      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.refFasta, prepareReferenceFiles.out.refFasta_primers, prepareReferenceFiles.out.bedfile, prepareReferenceFiles.out.primer_pairs, ch_blastndb, ch_krakendb)



}
