#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def help_message() {
    log.info"""
    ======================================================================================================================
    Nextflow workflow to get PofO assigned variants using long-read and Strand-seq.
    This nextflow pipeline is a comprehensive workflow to run all the steps 
    when long-read bam (with methylation tag) and Strand-seq fastqs are available.
    This workflow performs:
    1- Small variant calling using deepvariant or clair3 and phasing of them using whatshap or longphase.
    2- Large variant calling using sniffles.
    3- Adapter trimming, alignment, and QC of strand-seq data using cutadapt, 
        bowtie2, and ashleys-qc, respectively.
    4- Phasing small variants via strand-seq data using StrandPhaseR.
    5- Phasing and PofO assignment by also using LongPhase/WhatsHap phased blocks using patmat.py.
    ======================================================================================================================
    Usage:
    To run the nextflow: nextflow run patmat_workflow.nf -with-report -with-conda <specified optiones below>

    Mandatory arguments:
      --reference       Absolute path to reference genome. must be indexed by samtools faidx and bowtie2-build.
                        For bowtie2 indexing, you must use the reference name as base name 
                        (e.g. run bowtie2-build ref.fa ref.fa).
      --bam             Absolute Path to the long-read bam file with methylation tag.
      --output          Absolute Path to output directory
      --ashleys_model   Absolute Path to ashleys model pkl file.
      --strandseq_fq    Absolute path to the folder with strand-seq fastqs. NOTE: if your data are paired end, 
                        the paired reads must be in two different fastq files and the files must be named like 
                        UniqueID_R1.fastq and UniqueID_R2.fastq OR *_R1.fq and *_R2.fq (zipped .gz files are 
                        allowed. e.g. *_R1.fastq.gz). Unique ID will be a unique ID of each seqeunced cell.
    Optional arguments:
      --sample_id       Sample id. Will be also used as output prefix. Default is Sample
      --hifi            Select this option if long-read data is from PacBio
      --processes       Number of processes. Default is 10.
      --single          Select this if strand-seq is single-end and not paired
      --whatshap        Slelect this if you want to use whatshap instead of longphase
      --deepvar_model   <WGS|WES|PACBIO|ONT_R104|HYBRID_PACBIO_ILLUMINA>. Type of model to use 
                        for variant calling using deepvariant. Default is ONT_R104. 
      --clair3          Slelect this if you want to use clair3 instead of deepvariant
      --clair3_model    Specify the absolute path to clair3 model using this option if you want 
                        to use clair3 instead of deepvariant.
      --adapter_3R1     3 prime adapter of R1 for adapter trimming step. Default is
                        AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
      --adapter_5R1     5 prime adapter of R1 for adapter trimming step. Default is
                        AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT
      --adapter_3R2     3 prime adapter of R2 for adapter trimming step. Default is
                        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT
      --adapter_5R2     5 prime adapter of R2 for adapter trimming step. Default is
                        CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
      --help            Show help message
    """.stripIndent()
}

// Show help message
if (params.help){
    help_message()
    exit 0
}

// Parameters
params.reference = "reference file: not specified"
params.bam= "bam file: not specified"
params.output= "output directory: not specified"
params.ashleys_model= "ashleys model: not specified"
params.deepvar_model= "ONT_R104"
params.strandseq_fq= "strand-seq fastqs folder: not specified"
params.sample_id= "Sample"
params.processes= 10
params.adapter_3R1= "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
params.adapter_5R1= "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT"
params.adapter_3R2= "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"
params.adapter_5R2= "CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
params.clair3_model = "clair3 model: not specified"
params.clair3 = false
params.single = false
params.whatshap = false
params.hifi = false

def selected_params() {
    log.info"""
    --reference ${params.reference}
    --bam ${params.bam}
    --output ${params.output}
    --ashleys_model ${params.ashleys_model}
    --strandseq_fq ${params.strandseq_fq}
    --deepvar_model ${params.deepvar_model}
    --sample_id ${params.sample_id}
    --processes ${params.processes}
    --clair3 ${params.clair3}
    --clair3_model ${params.clair3_model}
    --single ${params.single}
    --whatshap ${params.whatshap}
    --hifi ${params.hifi}
    --adapter_3R1 ${params.adapter_3R1}
    --adapter_5R1 ${params.adapter_5R1}
    """.stripIndent()
}

println "\n\nParmateres For Running The Workflow (to see the help run the nextflow script with --help):"
if ( params.single ) {
    selected_params()
}
else {
    selected_params()
    println "--adapter_3R2 ${params.adapter_3R2}"
    println "--adapter_5R2 ${params.adapter_5R2}\n\n"
}

// Small Variant Calling
process small_variant_calling{
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/clair3_patmat-wf'
    publishDir "${params.output}/small_variant_results_${params.sample_id}", mode: 'copy'
    input:
        path(bam)
        path(bai)
        path(ref)
        path(fai)
    output:
        tuple path("*_passed_variants.vcf.gz"),
            path("{merge_output.vcf.gz,*_DeepVariant.vcf}")
    script:
        if ( params.clair3 ) {
            if ( params.hifi ) {
                """
                run_clair3.sh --bam_fn="$bam" \
                    --ref_fn="${ref}" --output="." \
                    --threads="${params.processes}" \
                    --platform="hifi" --model_path="${params.clair3_model}"
                gunzip -c merge_output.vcf.gz | awk '\$1~/^#/ || \$7=="PASS"' \
                | bgzip > "${params.sample_id}"_clair3_passed_variants.vcf.gz
                """
            }
            else {
                """
                    run_clair3.sh --bam_fn="$bam" \
                        --ref_fn="${ref}" --output="." \
                        --threads="${params.processes}" \
                        --platform="ont" --model_path="${params.clair3_model}"
                    gunzip -c merge_output.vcf.gz | awk '\$1~/^#/ || \$7=="PASS"' \
                    | bgzip > "${params.sample_id}"_clair3_passed_variants.vcf.gz
                    """
            }
        }
        else {
            """
            dir_b=\$(dirname "${ref}" | xargs realpath)
            dir_br=\$(realpath "${ref}" | rev | cut -d'/' -f2- | rev)
            dir_bb=\$(realpath "${bam}" | rev | cut -d'/' -f2- | rev)
            export SINGULARITY_CACHEDIR="\$dir_b"
            export SINGULARITY_TMPDIR="\$dir_b"
            BIN_VERSION="1.6.1"
            singularity pull docker://google/deepvariant:"\$BIN_VERSION"
            singularity run -B "\$dir_b":"\$dir_b" -B "\$dir_br":"\$dir_br" -B "\$dir_bb":"\$dir_bb" \
                docker://google/deepvariant:"\$BIN_VERSION" \
                /opt/deepvariant/bin/run_deepvariant \
                    --model_type="${params.deepvar_model}" \
                    --ref="${ref}" \
                    --reads="${bam}" \
                    --output_vcf="${params.sample_id}"_DeepVariant.vcf \
                    --num_shards="${params.processes}" \
                    --logging_dir="\$dir_b" \
                    --intermediate_results_dir="\$dir_b" \
                    --dry_run=false
            awk '\$1~/^#/ || \$7=="PASS"' "${params.sample_id}"_DeepVariant.vcf | \
            bgzip > "${params.sample_id}"_DeepVariant_passed_variants.vcf.gz
            """
        }
}



// SV Calling
process large_variant_calling{
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/patmat'
    publishDir "${params.output}/sniffles_results_${params.sample_id}", mode: 'copy'
    input:
        path(bam)
        path(bai)
        path(ref)
        path(fai)
        tuple path(clair3_pass_vcf), path(clair3_vcf)
    output:
        path("${params.sample_id}_sniffles_variants.vcf")
    script:
        """
        sniffles -i "${bam}" \
            --reference "${ref}" \
            -v "${params.sample_id}"_sniffles_variants.vcf
        """
}


// Phasing Long Reads
process phase_long_reads{
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/patmat'
    publishDir "${params.output}/LongReadPhased_${params.sample_id}", mode: 'copy'
    input:
        tuple path(pass_vcf), path(vcf)
        path(bam)
        path(bai)
        path(ref)
        path(fai)
        path(sniffles_vcf)
    output:
        tuple path("${params.sample_id}_*.vcf.gz"),
            path("${params.sample_id}_*.vcf.gz.tbi")
    script:
        if ( params.whatshap ) {
            """
            whatshap phase --ignore-read-groups \
                -o "${params.sample_id}"_WhatsHap.vcf \
                --reference="$ref" "$pass_vcf" "$bam"
            bgzip "${params.sample_id}"_WhatsHap.vcf
            tabix -p vcf "${params.sample_id}"_WhatsHap.vcf.gz
            """
        }
        else {
            if (params.hifi) {
                """
                longphase phase \
                    -s "$pass_vcf" -b "$bam" -r "$ref" --indels \
                    --pb -t "${params.processes}" -o "${params.sample_id}"_LongPhase
                bgzip "${params.sample_id}"_LongPhase.vcf
                tabix -p vcf "${params.sample_id}"_LongPhase.vcf.gz
                """
            }
            else {
                """
                longphase phase \
                    -s "$vcf" -b "$bam" -r "$ref" --indels \
                    --ont -t "${params.processes}" -o "${params.sample_id}"_LongPhase
                bgzip "${params.sample_id}"_LongPhase.vcf
                tabix -p vcf "${params.sample_id}"_LongPhase.vcf.gz
                """
            }
        }
}


// Adapter Trimming Of Strand-Seq Paired-End Reads
process strandseq_cutadapt_pair {
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/patmat'
    publishDir "${params.output}/cutadapt_${params.sample_id}", mode: 'copy'
    cpus "${params.processes}"
    input:
        tuple val(id),
            path(fastq_files)
        path(phased_long_reads_vcf)
    output:
        tuple val("${id}"),
            path("${id}_trimmed.1.fastq.gz"),
            path("${id}_trimmed.2.fastq.gz")
    script:
        def(read1, read2) = fastq_files
        """
        cutadapt \
            -a "${params.adapter_3R1}" \
            -g "${params.adapter_5R1}" \
            -A "${params.adapter_3R2}" \
            -G "${params.adapter_5R2}" \
            -o "${id}"_trimmed.1.fastq.gz \
            -p "${id}"_trimmed.2.fastq.gz \
            "${read1}" "${read2}" \
            -m 30 -q 15 -j "${params.processes}"
        """
}


// Adapter Trimming Of Strand-Seq Single-End Reads
process strandseq_cutadapt_single {
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/patmat'
    publishDir "${params.output}/cutadapt_${params.sample_id}", mode: 'copy'
    cpus "${params.processes}"
    input:
        path(fastq_file)
        path(phased_long_reads_vcf)
    output:
        path("*_trimmed.fastq.gz")
    script:
        """
        name=\$(echo "${fastq_file}" | rev | cut -d'/' -f1 | rev | sed 's/.fastq//g;s/.fq//g;s/.gz//g')
        cutadapt \
            -a "${params.adapter_3R1}" \
            -g "${params.adapter_5R1}" \
            -o "\$name"_trimmed.fastq.gz \
            "${fastq_file}" \
            -m 30 -q 15 -j "${params.processes}"
        """
}


// Alignment Of Strand-Seq Paired-End Reads
process strandseq_bowtie_pair {
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/patmat'
    publishDir "${params.output}/bowtie2_${params.sample_id}", mode: 'copy'
    cpus "${params.processes}"
    input:
        tuple val(id),
            path(fastq1),
            path(fastq2)
        path(bowtie_ref)
        path(bowtie_ref_index)
    output:
        path("${id}_MarkedDup.bam*")
    script:
        """
        chrs=\$(printf "chr%s " {1..22} X Y)
        bowtie2 --rg-id ${id} --rg "SM:${params.sample_id}" -x ${bowtie_ref} \
            -p "${params.processes}" -1 ${fastq1} -2 ${fastq2} | \
            samtools sort -@ "${params.processes}" -o "${id}"temp.bam
        samtools index -@ "${params.processes}" "${id}"temp.bam
        samtools view -h -F2052 -q10 "${id}"temp.bam \$chrs | \
            grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*_alt|@SQ.*chrM|@SQ.*chrEBV' | \
            samtools sort -@ "${params.processes}" -n | samtools fixmate -@ "${params.processes}" -O bam - - | \
            samtools view -bh -f1 | samtools sort -@ "${params.processes}" \
            -o "${id}".bam
        samtools index -@ "${params.processes}" "${id}".bam
        rm "${id}"temp.bam*
        picard MarkDuplicates -I "${id}".bam \
            -O "${id}"_MarkedDup.bam \
            -M "${id}"_MarkedDup.metrics
        rm "${id}".bam*
        samtools index -@ "${params.processes}" "${id}"_MarkedDup.bam
        """
}


// Alignment Of Strand-Seq Single-End Reads
process strandseq_bowtie_single {
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/patmat'
    publishDir "${params.output}/bowtie2_${params.sample_id}", mode: 'copy'
    cpus "${params.processes}"
    input:
        path(fastq)
        path(bowtie_ref)
        path(bowtie_ref_index)
    output:
        path("*_MarkedDup.bam*")
    script:
        """
        chrs=\$(printf "chr%s " {1..22} X Y)
        name=\$(echo "${fastq}" | rev | cut -d'/' -f1 | rev | sed 's/.fastq//g;s/.fq//g;s/.gz//g')
        bowtie2 --rg-id "\$name" --rg "SM:${params.sample_id}" -x ${bowtie_ref} -p "${params.processes}" -U ${fastq} | \
            samtools sort -@ "${params.processes}" -o "\$name"temp.bam
        samtools index -@ "${params.processes}" "\$name"temp.bam
        samtools view -h -F2052 -q10 "\$name"temp.bam \$chrs | \
            grep -v -E '@SQ.*chrUn|@SQ.*random|@SQ.*_alt|@SQ.*chrM|@SQ.*chrEBV' | \
            samtools sort -@ "${params.processes}" -o "\$name".bam
        samtools index -@ "${params.processes}" "\$name".bam
        rm "\$name"temp.bam* 
        picard MarkDuplicates -I "\$name".bam \
            -O "\$name"_MarkedDup.bam \
            -M "\$name"_MarkedDup.metrics
        rm "\$name".bam*
        samtools index -@ "${params.processes}" "\$name"_MarkedDup.bam
        """
}


// QC Strand-Seq Data
process ashley_qc{
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/ashleys_patmat-wf'
    publishDir "${params.output}/bowtie2_${params.sample_id}", mode: 'copy'
    input:
        path(strandseq_bam_folder)
        val(bowtie_out)
    output:
        tuple path("ashleys_pass"), path("ashleys_fail")
    script:
        """
        ashleys.py -j "${params.processes}" \
            features \
            -f ${strandseq_bam_folder} \
            -w 5000000 2000000 1000000 800000 600000 400000 200000 \
            -o ashleys_features.tsv

        ashleys.py -j "${params.processes}" \
            predict -p ashleys_features.tsv \
            -o ashleys_quality.txt \
            -m ${params.ashleys_model}

        mkdir ashleys_pass
        mkdir ashleys_fail

        awk 'NR>1 && \$2==1{ print \$1 }' ashleys_quality.txt | \
            xargs -i cp "${strandseq_bam_folder}"/{} "${strandseq_bam_folder}"/{}.bai ashleys_pass
        awk 'NR>1 && \$2==0{ print \$1 }' ashleys_quality.txt | \
            xargs -i cp "${strandseq_bam_folder}"/{} "${strandseq_bam_folder}"/{}.bai ashleys_fail

    """
}


// Phasing Of Strand-Seq Reads
process strandseq_phase {
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/patmat'
    publishDir "${params.output}/strandseq_phase_${params.sample_id}", mode: 'copy'
    input:
        tuple path(ashleys_pass_bam_folder),path(ashleys_fail_bam_folder)
        tuple path(clair3_pass_vcf), path(clair3_vcf)
    output:
        path("${params.sample_id}.strandseqphased.inv_aware.vcf")
    script:
        if (params.single){
            """
            strandseq_phase.R \
                -p FALSE \
                -i ${ashleys_pass_bam_folder} \
                -o . \
                -t "${params.processes}" \
                -n ${params.sample_id} \
                ${clair3_pass_vcf}
            """
        }
        else{
            """
            strandseq_phase.R \
                -p TRUE \
                -i ${ashleys_pass_bam_folder} \
                -o . \
                -t "${params.processes}" \
                -n ${params.sample_id} \
                ${clair3_vcf}
            """
        }
}


// Phasing Of Nanopore Reads and PofO Assignment
process patmat {
    tag "${params.sample_id}"
    conda '/projects/vakbari_prj/anaconda3/envs/patmat'
    publishDir "${params.output}/patmat_${params.sample_id}", mode: 'move'
    input:
        path(bam_file)
        path(bam_index)
        tuple path(longphase_vcf), path(longphase_index)
        path(sniffles_vcf)
        path(strandseq_phased_vcf)
        path(reference_genome)
        path(reference_index)
    output:
        path("${params.sample_id}*")
    script:
        if (params.hifi){
            """
            patmat.py \
                -v ${longphase_vcf} -ph \
                -b ${bam_file} -stv ${strandseq_phased_vcf} \
                -o "${params.sample_id}" \
                -p ${params.processes} -sv_vcf ${sniffles_vcf} \
                -ref ${reference_genome} -pb -is
            """
        }
        else{
            """
            patmat.py \
                -v ${longphase_vcf} -ph \
                -b ${bam_file} -stv ${strandseq_phased_vcf} \
                -o "${params.sample_id}" \
                -p ${params.processes} -sv_vcf ${sniffles_vcf} \
                -ref ${reference_genome} -is
            """
        }
}

workflow {
    Channel
        .fromFilePairs("${params.strandseq_fq}/*_{R1,R2}.{fastq,fq,fastq.gz,fq.gz}")
        .set{strandseq_fqs_pair}
    Channel
        .fromPath("${params.strandseq_fq}/*")
        .set{strandseq_fqs_single}
    Channel
        .fromPath("${params.reference}.*").collect()
        .set{ref_indexes}
    small_variant_calling(params.bam,file("${params.bam}.{bai,crai}"),params.reference,
                            "${params.reference}.fai")
    large_variant_calling(params.bam,file("${params.bam}.{bai,crai}"),params.reference,
                            "${params.reference}.fai",
                            small_variant_calling.out)
    phase_long_reads(small_variant_calling.out,params.bam,file("${params.bam}.{bai,crai}"),
                            params.reference,"${params.reference}.fai", 
                            large_variant_calling.out)
    if ( params.single ) {
        strandseq_cutadapt_single(strandseq_fqs_single,
                                    phase_long_reads.out)
        strandseq_bowtie_single(strandseq_cutadapt_single.out,
                                params.reference,ref_indexes)
        ashley_qc("${params.output}/bowtie2_${params.sample_id}",
                    strandseq_bowtie_single.out.collect(flat:true))
    }
    else {
        strandseq_cutadapt_pair(strandseq_fqs_pair,
                                phase_long_reads.out)
        strandseq_bowtie_pair(strandseq_cutadapt_pair.out,
                            params.reference, ref_indexes)
        ashley_qc("${params.output}/bowtie2_${params.sample_id}",
                    strandseq_bowtie_pair.out.collect(flat:true))
    }
    strandseq_phase(ashley_qc.out,
                    small_variant_calling.out)
    patmat(params.bam, file("${params.bam}.{bai,crai}"), phase_long_reads.out,
            large_variant_calling.out, strandseq_phase.out, 
            params.reference,"${params.reference}.fai")
}


