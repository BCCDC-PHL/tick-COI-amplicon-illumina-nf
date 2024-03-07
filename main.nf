#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// include modules
include {printHelp} from './modules/help.nf'

// import subworkflows
include {coiIllumina} from './workflows/illuminaCOI.nf'

if (params.help){
    printHelp()
    exit 0
}

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

if ( !params.directory ) {
  println("Please supply a directory containing fastqs or CRAMs with --directory.")
  println("Use --help to print help")
  System.exit(1)
}

if ( (params.bed && ! params.ref) || (!params.bed && params.ref) ) {
  println("--bed and --ref must be supplied together")
  System.exit(1)
}

if ( !params.ref_with_primers ) {
  println("Please supply a reference fasta file with primers")
  println("Use --help to print help")
  System.exit(1)
}

if ( ! params.prefix ) {
     println("Please supply a prefix for your output files with --prefix")
     println("Use --help to print help")
     System.exit(1)
} else {
     if ( params.prefix =~ /\// ){
         println("The --prefix that you supplied contains a \"/\", please replace it with another character")
         System.exit(1)
     }
} 


// main workflow
workflow {

  Channel.fromFilePairs( params.fastqSearchPath, flat: true)
	          .filter{ !( it[0] =~ /Undetermined/ ) }
	          .set{ ch_filePairs }


  Channel.fromPath(new File(params.db).getAbsolutePath())
		  .set{ ch_blastndb }

  Channel.fromPath(new File(params.kraken_db).getAbsolutePath())
                  .set{ ch_krakendb }

  main:
   coiIllumina(ch_filePairs, ch_blastndb, ch_krakendb)

}
