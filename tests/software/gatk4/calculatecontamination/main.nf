#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CALCULATECONTAMINATION } from '../../../../software/gatk4/calculatecontamination/main.nf' addParams( options: [:] )

workflow test_gatk4_calculatecontamination {
    
    def input = []
    input = [ [ id:'test'], // meta map
              [file("${launchDir}/tests/data/gatk/7_tumor_getpileupsummaries.table", checkIfExists: true) ]

    GATK4_CALCULATECONTAMINATION ( input )
}
