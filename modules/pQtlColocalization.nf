#!/bin/bash nextflow


process ExtractGenesForPQtlAnalysis {
    scratch true

    input:
        tuple val(id), val(ensembl), path(pqtl)
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
        tuple val(id), val(ensembl), path(pqtl)

    output:
        tuple val(id), path("concatenated.processed.txt")

    shell:
        variant_id = "ID"
        chr = "CHROM"
        pos = "GENPOS"
        effect_allele = "ALLELE1"
        other_allele = "ALLELE0"
        effect = "BETA"
        standard_error = "SE"
        mlog10_p_value = "LOG10P"
        allele_frequency = "A1FREQ"
        liftover = "-"
        '''
        mkdir untar
        tar -C untar -xvf !{pqtl}

        files=(untar/*/*.gz)
        total=${#files[@]}
        i=0

        for f in "${files[@]}"; do
            i=$(( i + 1 ))
            AdjustGwasFile.R -i ${f} \
            -id !{variant_id} -chr !{chr} -bp !{pos} -ea !{effect_allele} -oa !{other_allele} \
            -b !{effect} -se !{standard_error} -mlp !{mlog10_p_value} -af !{allele_frequency} -l !{liftover} \
            -out !{id}_${i}_processed.txt

        done

        head -n 1 "!{id}_0_processed.txt" > "concatenated.processed.txt"
        tail -n +2 "!{id}_*_processed.txt" >> "concatenated.processed.txt"

        '''
}

process ComparePqtlAndEqtl {
    input:
        tuple val(id), path(pqtl), path(eqtl)

    output:
        path "output"

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

    main:
        extracted_ch = ExtractGenesForPQtlAnalysis(
            pqtl_ch, eqtl_ch)

        joined_preprocessed_ch = AdjustPQtlFile(pqtl_ch).join(extracted_ch)

        comparison_output_ch = ComparePqtlAndEqtl(joined_preprocessed_ch)

    emit:
        comparison_output_ch
}