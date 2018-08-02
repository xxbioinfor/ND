package org.synergylab.neoantigenDiscovery.utils

//设定文件基础路径、参考文件路径、软件路径等信息

object ProjectPath{

    val baseDir = "/localdisk/"
    val softwareDir = baseDir + "software/"
    val referenceDir = baseDir + "reference/"
    val msDir = "E:/MS" //质谱位于wins系统
    val jobDir = baseDir + "jobs/"
    val thread = 20

    //userConn
    val linuxIP = "192.168.6.203"
    val linuxName = "root"
    val linuxPwd = "rlibs402"
    val winIP = "192.168.6.204"
    val winName = "rlibs"
    val winPwd = "rlibs402"

    //ms

    //jobs
    val cancerType = "lung"
    val cancerSampleName_1 = "180313002ML_S3_L008_R1_001.fastq.gz"
    val cancerSampleName_2 = "180313002ML_S3_L008_R2_001.fastq.gz"
    val controlSampleName_1 = "180313001ML_S14_L008_R1_001.fastq.gz"
    val controlSampleName_2 = "180313001ML_S14_L008_R2_001.fastq.gz"
    val sampleID = getSubUtilSimple(cancerSampleName_1,"(.*?)_.*")
    val sampleBaseDir = jobDir + "lung/" + sampleID + "/"
    val cancerSampleDir = sampleBaseDir + "cancer/"
    val controlSampleDir = sampleBaseDir + "control/"
    val rawCancerSampleDir = cancerSampleDir + "rawdata/"
    val cleanCancerSampleDir = cancerSampleDir + "cleandata/"
    val rawControlSampleDir = controlSampleDir + "rawdata/"
    val cleanControlSampleDir = controlSampleDir + "cleandata/"
    //val cancerSampleFilePrefix = cancerSampleDir + sampleID
    //val controlSampleFilePrefix = controlSampleDir + sampleID
    val rawCancerSequence_1 = rawCancerSampleDir + cancerSampleName_1
    val rawCancerSequence_2 = rawCancerSampleDir + cancerSampleName_2
    val rawControlSequence_1 = rawControlSampleDir + controlSampleName_1
    val rawControlSequence_2 = rawControlSampleDir + controlSampleName_2

    //创建辅助jobs路径

    //software
    val fastqc = softwareDir + "FastQC/fastqc"
    val trimmomatic = softwareDir + "Trimmomatic-0.36/trimmomatic-0.36.jar"
    val hisat2 = softwareDir + "hisat2-2.1.0/hisat2"
    val samtools = softwareDir + "samtools/samtools"
    val trinity = softwareDir + "trinityrnaseq-Trinity-v2.5.1/"
    val picard = softwareDir + "picard.jar"
    val GATK = softwareDir + "GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar"
    val vep = softwareDir + "ensembl-vep/vep"
    val RSEM = softwareDir + "RSEM-1.3.0/rsem-calculate-expression"
    val bowtie2 = softwareDir + "bowtie2-2.3.3.1-linux-x86_64/bowtie2"
    val seq2HLA = softwareDir + "seq2HLA2.2/seq2HLA.py"
    val netMHC = softwareDir + "netMHC-4.0/netMHC"
    val gmap = softwareDir + "gmap-2018-03-01/"
    //质谱工具未列出


    //reference
    val hisat2Index = referenceDir + "hisat2_index/hisat2_index"
    val gtfReference = referenceDir + "Homo_sapiens.GRCh38.91.gtf" //有三个版本，Chr、chr、空，后续修改
    val hg38Reference = referenceDir + "Homo_sapiens_assembly38.fasta"
    val faReference = referenceDir + "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    val vcfReference = referenceDir + "known_vcf/"
    val mills_1000G = vcfReference + "Mills_and_1000G_gold_standard.indels.hg38.vcf"
    val omni_1000G = vcfReference + "1000G_omni2.5.hg38.vcf"
    val dbsnp = vcfReference + "dbsnp_146.hg38.vcf"
    val hapmap = vcfReference + "hapmap_3.3.hg38.vcf"
    val cosmic = vcfReference + "CosmicCodingMuts.vcf"  //可能用不到
    val adapter = referenceDir + "adapter.fasta"
    val vepReference = referenceDir + ".vep"

}

