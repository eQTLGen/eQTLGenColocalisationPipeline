#!/bin/bash nextflow


process ExtractGenesForPQtlAnalysis {
    scratch true

    input:
        tuple val(id), path(pqtl), val(ensembl)
        path input

    output:
        tuple val(id), path("extracted*")

    shell:
        '''
        mkdir tmp_eqtls

        cp -r "!{input}/phenotype=!{ensembl}" tmp_eqtls/

        extract_parquet_results.py \
            --input-file tmp_eqtls \
            --genes !{ensembl} \
            --cols "+p_value" \
            --output-prefix extracted

        rm -r tmp_eqtls
        '''
}

process AdjustPQtlFile {

    input:
        tuple val(id), path(pqtl), val(ensembl)

    output:
        tuple val(id), path("*.processed.txt").collectFile(skip:1, keepHeader:true)

    shell:
        def variant_id = 3
        def chr = 1
        def pos = 2
        def effect_allele = 5
        def other_allele = 4
        def effect = 10
        def standard_error = 11
        def p_value = 13
        def allele_frequency = 6
        def liftover = ""
        '''
        tar -C untar -xvf !{pqtl}

        files=(untar/*/*.gz)
        total=${#files[@]}
        i=0

        for f in "${files[@]}"; do
            i=$(( i + 1 ))
            Rscript --vanilla $baseDir/bin/AdjustGwasFile.R !{gwas} \
            !{variant_id} !{chr} !{pos} !{effect_allele} !{other_allele} !{effect} !{standard_error} !{p_value} !{allele_frequency} !{liftover} \
            !{id}_${i}_processed.txt

        done
        '''

}

process ComparePqtlAndEqtl {
    input:
        tuple val(id), path(pqtl), path(eqtl)

    output:

    script:
        """
        # Make manhattan plot for both
        # For each
        comparePqtlToEqtl.R --id ${id} --protein ${pqtl} --expression ${eqtl}
        """
}

workflow PQTL_COMPARISON {
    take:
        pqtl_ch
        eqtl_ch
        posterior_threshold
        cs_threshold
        output_cs_pip

    main:
        extracted_ch = ExtractGenesForPQtlAnalysis(
            pqtl_ch, eqtl_ch)

        joined_preprocessed_ch = AdjustPQtlFile(pqtl_ch).join(extracted_ch)

        ComparePqtlAndEqtl(joined_preprocessed_ch)

    emit:
        cs = hypr_coloc_results.cs
        pips = hypr_coloc_results.pips
}