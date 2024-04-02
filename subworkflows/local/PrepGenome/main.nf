include { BWAMEM2_INDEX }                       from '../../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX }                      from '../../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY }     from '../../../modules/nf-core/picard/createsequencedictionary/main'


workflow PrepGenome {
    take:
    fasta       // channel: [mandatory] [fasta]
    fasta_fai

    main:
    // Initialize output channels
    versions = Channel.empty()

    // Map fasta to [meta, fasta] channel
    fasta = fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] }

    // Only runs if bwa2_index is not provided in params
    BWAMEM2_INDEX(fasta)

    SAMTOOLS_FAIDX(fasta, [[id:'null'], []])
    PICARD_CREATESEQUENCEDICTIONARY(fasta)

    versions = versions.mix(BWAMEM2_INDEX.out.versions)

    emit:
    bwa_index = BWAMEM2_INDEX.out.index.map{ meta, index -> [index] }.collect()     // path: bwamem2/*
    fasta_fai = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }
    fasta_dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict.map{ meta, dict -> [dict] }
    versions                                                                        // channel: [ versions.yml ]
}

