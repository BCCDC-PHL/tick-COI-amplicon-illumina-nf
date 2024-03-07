def printHelp() {
  log.info"""
  Usage:
    nextflow run BCCDC-PHL/tick-COI-amplicon-illumina-nf -profile conda --cached ~/.conda/envs --prefix [prefix] [workflow-options]

  Description:
    Turn Illumina sequencing reads into consensus sequences 

  Nextflow arguments (single DASH):
    -profile                      Allowed values: conda
 
  Workflow options:
    Mandatory:
      --prefix                    A (unique) string prefix for output files.
                                  Sequencing run name is a good choice e.g DDMMYY_MACHINEID_RUN_FLOWCELLID.
      --directory                 Path to a directory containing paired-end Illumina reads. 
                                  Reads will be found and paired RECURSIVELY beneath this directory.
      --db			  Path to blastn database (e.g., in-house tick COI database or supply your specific database)
      --kraken_db                 Path to kraken2 database
    Optional:
      --outdir                    Output directory (Default: ./results) 
      --bed                       Path to primer bed file, also requires --ref-with-primers
      --ref                       Path to reference fasta file (without primers)
      --ref_with_primers          Path to reference fasta file (with primers), also requires --bed
      --primer_pairs_tsv          File showing which primers are paired.
      --composite_ref             Human_and_composite_ref sequence
      --length                    Threshold for keeping contigs from denovoAssembly (default:600)
      --qualThreshold             Sliding window quality threshold for keeping reads (default: 20)
      --cache                     Specify that Conda environments should be cached in a specific directory (example: ~/.conda/envs) 

  """.stripIndent()
}
