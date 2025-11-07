# 环境配置方面有一点点小bug，需要调一下；
# 初步发现的需要隔离环境的是sra-tools，我用的服务器不隔离会下载旧版本。

set -e 
set -o pipefail 
echo "GATK_ana: "
echo "start time: $(date)"

PROJECT_DIR="/root/autodl-tmp/my_first_project"
DATA_DIR="${PROJECT_DIR}/data"
REF_DIR="${PROJECT_DIR}/reference"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs" 
GENOMICSDB_WORKSPACE="${RESULTS_DIR}/genomicsdb_workspace"

REF_GENOME="hg38.fa" 
SAMPLES=("final_test") 

DBSNP="dbsnp_146.hg38.vcf.gz"
MILLS_INDELS="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_INDELS="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

: <<'VQSR'
# 多样本分析取消禁用
HAPMAP="hapmap_3.3.hg38.vcf.gz"
OMNI="1000G_omni2.5.hg38.vcf.gz"
PHASE1_SNPS="1000G_phase1.snps.high_confidence.hg38.vcf.gz"
VQSR

THREADS=16 
MEM_GB=100 

# log
mkdir -p "${RESULTS_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${GENOMICSDB_WORKSPACE}"
# 索引
echo "indexing: "
if [ ! -f "${REF_DIR}/${REF_GENOME}.fai" ]; then
    echo "process1: Samtools"
    samtools faidx "${REF_DIR}/${REF_GENOME}"
    gatk CreateSequenceDictionary -R "${REF_DIR}/${REF_GENOME}"
else
    echo "skip_samtools: existed index"
fi

if [ ! -f "${REF_DIR}/${REF_GENOME}.bwt" ]; then
    echo "process2: BWA"
    bwa index "${REF_DIR}/${REF_GENOME}"
else
    echo "skip_BWA: existed index"
fi
echo "indexing done:) time: $(date)"


# 质控
echo "PROCESSING: QUALITY CONTROL"
GVCF_FILES=() 

for SAMPLE_ID in "${SAMPLES[@]}"; do
    echo "processing: ${SAMPLE_ID}"

    FASTQ_R1="${DATA_DIR}/${SAMPLE_ID}_R1.fastq.gz"
    FASTQ_R2="${DATA_DIR}/${SAMPLE_ID}_R2.fastq.gz"
    RAW_BAM="${RESULTS_DIR}/${SAMPLE_ID}.raw.bam"
    SORTED_BAM="${RESULTS_DIR}/${SAMPLE_ID}.sorted.bam"
    DEDUP_BAM="${RESULTS_DIR}/${SAMPLE_ID}.dedup.bam"
    RECAL_TABLE="${RESULTS_DIR}/${SAMPLE_ID}.recal_table"
    ANALYSIS_READY_BAM="${RESULTS_DIR}/${SAMPLE_ID}.analysis_ready.bam"
    GVCF_FILE="${RESULTS_DIR}/${SAMPLE_ID}.g.vcf.gz"
    
    GVCF_FILES+=("$GVCF_FILE")

    # read_group
    RG_ID="group_${SAMPLE_ID}"
    RG_PU="unit1"
    RG_SM="${SAMPLE_ID}"
    RG_PL="ILLUMINA"
    RG_LB="lib1"

    # BWA
    echo "processing: BWA_alignment"
    bwa mem -t ${THREADS} -R "@RG\tID:${RG_ID}\tPU:${RG_PU}\tSM:${RG_SM}\tPL:${RG_PL}\tLB:${RG_LB}" \
        "${REF_DIR}/${REF_GENOME}" "${FASTQ_R1}" "${FASTQ_R2}" | \
        samtools sort -@ ${THREADS} -o "${SORTED_BAM}" -
    # PCR
    echo "processing: mark_duplicates"
    gatk --java-options "-Xmx${MEM_GB}G" MarkDuplicates \
        -I "${SORTED_BAM}" \
        -O "${DEDUP_BAM}" \
        -M "${RESULTS_DIR}/${SAMPLE_ID}.dedup_metrics.txt"
    # BQSR
    echo "processing: BQSR_recalibration"
    gatk --java-options "-Xmx${MEM_GB}G" BaseRecalibrator \
        -R "${REF_DIR}/${REF_GENOME}" \
        -I "${DEDUP_BAM}" \
        --known-sites "${REF_DIR}/${DBSNP}" \
        --known-sites "${REF_DIR}/${MILLS_INDELS}" \
        --known-sites "${REF_DIR}/${KNOWN_INDELS}" \
        -O "${RECAL_TABLE}"

    echo "applying_BQSR"
    gatk --java-options "-Xmx${MEM_GB}G" ApplyBQSR \
        -R "${REF_DIR}/${REF_GENOME}" \
        -I "${DEDUP_BAM}" \
        -bqsr "${RECAL_TABLE}" \
        -O "${ANALYSIS_READY_BAM}"
    # HaplotypeCaller
    echo "HaplotypeCaller"
    gatk --java-options "-Xmx${MEM_GB}G" HaplotypeCaller \
        -R "${REF_DIR}/${REF_GENOME}" \
        -I "${ANALYSIS_READY_BAM}" \
        -O "${GVCF_FILE}" \
        -ERC GVCF
        
    echo "done: ${SAMPLE_ID}"
done
echo "QUALITY CONTROL DONE"

# 合并&分型
echo "PROCESSING: MERGE AND GENOTYPE"
# V参
GVCF_ARGS=""
for gvcf in "${GVCF_FILES[@]}"; do
    GVCF_ARGS+=" -V ${gvcf}"
done

# wspace clean
echo "cleaning workspace"
rm -rf "${GENOMICSDB_WORKSPACE}"

gatk --java-options "-Xmx${MEM_GB}G" GenomicsDBImport \
    -L "${REF_DIR}/hg38.intervals.list" \
    --genomicsdb-workspace-path "${GENOMICSDB_WORKSPACE}" \
    --batch-size 50 \
    --reader-threads ${THREADS} \
    ${GVCF_ARGS}

# GenotypeGVCFs
echo "GenotypeGVCFs: "
RAW_VCF="${RESULTS_DIR}/cohort.raw.vcf.gz"
gatk --java-options "-Xmx${MEM_GB}G" GenotypeGVCFs \
    -R "${REF_DIR}/${REF_GENOME}" \
    -V "gendb://${GENOMICSDB_WORKSPACE}" \
    -O "${RAW_VCF}"
echo "done"


echo "HARD FILTERING"
FILTERED_VCF="${RESULTS_DIR}/cohort.filtered.vcf.gz"

echo "seperate SNP"
RAW_SNPS_VCF="${RESULTS_DIR}/cohort.raw.snps.vcf.gz"
gatk --java-options "-Xmx${MEM_GB}G" SelectVariants \
    -V "${RAW_VCF}" \
    -O "${RAW_SNPS_VCF}" \
    -select-type SNP
echo "seperate INDEL"
RAW_INDELS_VCF="${RESULTS_DIR}/cohort.raw.indels.vcf.gz"
gatk --java-options "-Xmx${MEM_GB}G" SelectVariants \
    -V "${RAW_VCF}" \
    -O "${RAW_INDELS_VCF}" \
    -select-type INDEL

echo "SNP filtering"
FILTERED_SNPS_VCF="${RESULTS_DIR}/cohort.filtered.snps.vcf.gz"
gatk --java-options "-Xmx${MEM_GB}G" VariantFiltration \
    -V "${RAW_SNPS_VCF}" \
    -O "${FILTERED_SNPS_VCF}" \
    -filter "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
    -filter-name "gatk_snp_filter"

echo "indel filtering"
FILTERED_INDELS_VCF="${RESULTS_DIR}/cohort.filtered.indels.vcf.gz"
gatk --java-options "-Xmx${MEM_GB}G" VariantFiltration \
    -V "${RAW_INDELS_VCF}" \
    -O "${FILTERED_INDELS_VCF}" \
    -filter "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
    -filter-name "gatk_indel_filter"
# hebing
echo "merge"
gatk --java-options "-Xmx${MEM_GB}G" MergeVcfs \
    -I "${FILTERED_SNPS_VCF}" \
    -I "${FILTERED_INDELS_VCF}" \
    -O "${FILTERED_VCF}"

echo "fin"

echo "${FILTERED_VCF}"
echo "fin_t: $(date)"

# 多样本VQSR，取消禁用
: <<'VQSRana'
RECAL_SNP_FILE="${RESULTS_DIR}/cohort.snps.recal"
TRANCHES_SNP_FILE="${RESULTS_DIR}/cohort.snps.tranches"
RECAL_INDEL_FILE="${RESULTS_DIR}/cohort.indels.recal"
TRANCHES_INDEL_FILE="${RESULTS_DIR}/cohort.indels.tranches"
FILTERED_VCF="${RESULTS_DIR}/cohort.filtered.vcf.gz"

echo "building SNP model"
gatk --java-options "-Xmx${MEM_GB}G" VariantRecalibrator \
    -V "${RAW_VCF}" \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${REF_DIR}/${HAPMAP}" \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 "${REF_DIR}/${OMNI}" \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 "${REF_DIR}/${PHASE1_SNPS}" \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${REF_DIR}/${DBSNP}" \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O "${RECAL_SNP_FILE}" \
    --tranches-file "${TRANCHES_SNP_FILE}"

echo "building indel model"
gatk --java-options "-Xmx${MEM_GB}G" VariantRecalibrator \
    -V "${RAW_VCF}" \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 "${REF_DIR}/${MILLS_INDELS}" \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${REF_DIR}/${DBSNP}" \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    -O "${RECAL_INDEL_FILE}" \
    --tranches-file "${TRANCHES_INDEL_FILE}"

echo "applying VQSR_filter"
TEMP_VCF1="${RESULTS_DIR}/cohort.temp1.vcf.gz"
gatk --java-options "-Xmx${MEM_GB}G" ApplyVQSR \
    -V "${RAW_VCF}" \
    -O "${TEMP_VCF1}" \
    --recal-file "${RECAL_INDEL_FILE}" \
    --tranches-file "${TRANCHES_INDEL_FILE}" \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode INDEL

gatk --java-options "-Xmx${MEM_GB}G" ApplyVQSR \
    -V "${TEMP_VCF1}" \
    -O "${FILTERED_VCF}" \
    --recal-file "${RECAL_SNP_FILE}" \
    --tranches-file "${TRANCHES_SNP_FILE}" \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode SNP
echo "done"

VQSRana
