package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//由bam文件，经由Picard、GATK、MuTect2、VEP确认非同义体细胞突变

class NSSomaticMutation(){
    //parameters

    //source data
    val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)

    //class properties
    //val somaticMutation ＝ generateSomaticMutation()

    //function
    fun picardCommand(bamFile: String,sampleName: String,sampleID: String,outDir: String): String {
       rec.execute("java -Xms20g -jar "+ProjectPath.picard+" AddOrReplaceReadGroups INPUT="+bamFile+" OUTPUT="+
               outDir+sampleID+"_sorted.bam SORT_ORDER=coordinate RGPL=illumina RGLB="+sampleName+" RGPU=temp RGSM="+sampleName)
       rec.execute("java -Xms20g -jar "+ProjectPath.picard+" MarkDuplicates I="+outDir+sampleID+"_sorted.bam O="
               +outDir+sampleID+"_duplicateRemoved.bam M=" +outDir+sampleID+"_marked_dup_metrics.txt REMOVE_DUPLICATES=true AS=true")
       rec.execute("java -Xms20g -jar "+ProjectPath.picard+" ReorderSam I="+outDir+sampleID +"_duplicateRemoved.bam O="
               +outDir+sampleID+"_Reorder.bam REFERENCE=" +ProjectPath.hg38Reference)
       rec.execute("java -Xms20g -jar "+ProjectPath.picard+" BuildBamIndex INPUT="+outDir+sampleID+"_Reorder.bam")

       return outDir+sampleID+"_Reorder.bam"
    }

    fun gatkCommand(reorderBamFile: String,sampleID: String,outDir: String): String {
       rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T SplitNCigarReads -R "+ProjectPath.hg38Reference+" -I "
               +reorderBamFile+" -o "+outDir+sampleID+"_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS")
       rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T RealignerTargetCreator -R "+ProjectPath.hg38Reference
               +" -I "+outDir+sampleID+"_split.bam --known "+ProjectPath.omni_1000G +" --known "+ProjectPath.dbsnp
               +" --known "+ProjectPath.hapmap +" --filter_reads_with_N_cigar -o " +outDir+sampleID+".intervals")
       rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T IndelRealigner -R "+ProjectPath.hg38Reference+" -I "
               +outDir+sampleID+"_Reorder.bam -targetIntervals " +outDir+sampleID+".intervals -known "+ProjectPath.omni_1000G
               +" -known "+ProjectPath.dbsnp+" -known "+ProjectPath.hapmap +" --filter_reads_with_N_cigar -o "
               +outDir+sampleID+"_IndelRealign.bam")
       rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T BaseRecalibrator -R "+ProjectPath.hg38Reference+" -I "
               +outDir+sampleID+"_IndelRealign.bam -knownSites "+ProjectPath.hapmap +" -o "+outDir+sampleID+"_Recal.table -nct "+ProjectPath.thread)
       rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T PrintReads -R "+ProjectPath.hg38Reference+" -I "
               +outDir+sampleID+"_IndelRealign.bam -BQSR "+outDir +sampleID+"_Recal.table -o "+outDir+sampleID
               +"_Recal.bam -nct "+ProjectPath.thread)

       return outDir+sampleID+"_Recal.bam"
    }

    fun mutect2Command(cancerRecalBamFile: String,controlRecalBamFile: String,outDir: String): String {
       rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T MuTect2 -R "+ProjectPath.hg38Reference+" -I:tumor "
               +cancerRecalBamFile+" -I:normal "+controlRecalBamFile+" -o " +outDir
               +"somatic_mutations.vcf -nct "+ProjectPath.thread)
       return outDir+"somatic_mutations.vcf"
    }

    fun vepCommand(mutationFile: String,outDir: String): String {
       rec.execute(ProjectPath.vep+" --cache --dir "+ProjectPath.vepReference+" --format vcf --input_file "
               +mutationFile+" --output_file "+outDir+"somatic_mutations_pass_vep.vcf")
        return outDir+"somatic_mutations_pass_vep.vcf"
    }
   /*
    //function
     fun generateSomaticMutation(){
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)

        //picard
        //cancer
        rec.execute("java -Xms20g -jar "+ProjectPath.picard+" AddOrReplaceReadGroups INPUT="+ProjectPath.cancerSampleDir
                +ProjectPath.sampleID+"_cancer.bam OUTPUT="+ProjectPath.cancerSampleDir+ProjectPath.sampleID
                +"_cancer_sorted.bam SORT_ORDER=coordinate RGPL=illumina RGLB=cancer RGPU=temp RGSM=cancer")
        rec.execute("java -Xms20g -jar "+ProjectPath.picard+" MarkDuplicates I="+ProjectPath.cancerSampleDir
                +ProjectPath.sampleID+"_cancer_sorted.bam O="+ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_duplicateRemoved.bam M="
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_marked_dup_metrics.txt REMOVE_DUPLICATES=true AS=true")
        rec.execute("java -Xms20g -jar "+ProjectPath.picard+" ReorderSam I="+ProjectPath.cancerSampleDir+ProjectPath.sampleID
                +"_cancer_duplicateRemoved.bam O="+ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_Reorder.bam REFERENCE="
                +ProjectPath.hg38Reference)
        rec.execute("java -Xms20g -jar "+ProjectPath.picard+" BuildBamIndex INPUT="+ProjectPath.cancerSampleDir
                +ProjectPath.sampleID+"_cancer_Reorder.bam")
        //control
        rec.execute("java -Xms20g -jar "+ProjectPath.picard+" AddOrReplaceReadGroups INPUT="+ProjectPath.controlSampleDir
                +ProjectPath.sampleID+"_control.bam OUTPUT="+ProjectPath.controlSampleDir+ProjectPath.sampleID
                +"_control_sorted.bam SORT_ORDER=coordinate RGPL=illumina RGLB=control RGPU=temp RGSM=control")
        rec.execute("java -Xms20g -jar "+ProjectPath.picard+" MarkDuplicates I="+ProjectPath.controlSampleDir
                +ProjectPath.sampleID+"_control_sorted.bam O="+ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_duplicateRemoved.bam M="
                +ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_marked_dup_metrics.txt REMOVE_DUPLICATES=true AS=true")
        rec.execute("java -Xms20g -jar "+ProjectPath.picard+" ReorderSam I="+ProjectPath.controlSampleDir+ProjectPath.sampleID
                +"_control_duplicateRemoved.bam O="+ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_Reorder.bam REFERENCE="
                +ProjectPath.hg38Reference)
        rec.execute("java -Xms20g -jar "+ProjectPath.picard+" BuildBamIndex INPUT="+ProjectPath.controlSampleDir
                +ProjectPath.sampleID+"_control_Reorder.bam")

        //GATK
        //cancer
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T SplitNCigarReads -R "+ProjectPath.hg38Reference+" -I "
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_Reorder.bam -o "+ProjectPath.cancerSampleDir
                +ProjectPath.sampleID+"_cancer_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS")
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T RealignerTargetCreator -R "+ProjectPath.hg38Reference
                +" -I "+ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_split.bam --known "+ProjectPath.omni_1000G
                +" --known "+ProjectPath.dbsnp+" --known "+ProjectPath.hapmap +" --filter_reads_with_N_cigar -o"
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer.intervals")
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T IndelRealigner -R "+ProjectPath.hg38Reference+" -I "
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_Reorder.bam -targetIntervals "
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer.intervals -known "+ProjectPath.omni_1000G
                +" -known "+ProjectPath.dbsnp+" -known "+ProjectPath.hapmap +" --filter_reads_with_N_cigar -o "
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_IndelRealign.bam")
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T BaseRecalibrator -R "+ProjectPath.hg38Reference+" -I "
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_IndelRealign.bam -knownSites "+ProjectPath.hapmap
                +" -o "+ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_Recal.table -nct "+ProjectPath.thread)
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T PrintReads -R "+ProjectPath.hg38Reference+" -I "
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_IndelRealign.bam -BQSR "+ProjectPath.cancerSampleDir
                +ProjectPath.sampleID+"_cancer_Recal.table -o "+ProjectPath.cancerSampleDir+ProjectPath.sampleID
                +"_cancer_Recal.bam -nct "+ProjectPath.thread)

        //control
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T SplitNCigarReads -R "+ProjectPath.hg38Reference+" -I "
                +ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_Reorder.bam -o "+ProjectPath.controlSampleDir
                +ProjectPath.sampleID+"_control_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS")
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T RealignerTargetCreator -R "+ProjectPath.hg38Reference
                +" -I "+ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_split.bam --known "+ProjectPath.omni_1000G
                +" --known "+ProjectPath.dbsnp+" --known "+ProjectPath.hapmap+" --filter_reads_with_N_cigar -o"
                +ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control.intervals")
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T IndelRealigner -R "+ProjectPath.hg38Reference+" -I "
                +ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_Reorder.bam -targetIntervals "
                +ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control.intervals -known "+ProjectPath.omni_1000G
                +" -known "+ProjectPath.dbsnp+" -known "+ProjectPath.hapmap+" --filter_reads_with_N_cigar -o "
                +ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_IndelRealign.bam")
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T BaseRecalibrator -R "+ProjectPath.hg38Reference+" -I "
                +ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_IndelRealign.bam -knownSites "+ProjectPath.hapmap
                +" -o "+ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_Recal.table -nct "+ProjectPath.thread)
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T PrintReads -R "+ProjectPath.hg38Reference+" -I "
                +ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control_IndelRealign.bam -BQSR "+ProjectPath.controlSampleDir
                +ProjectPath.sampleID+"_control_Recal.table -o "+ProjectPath.controlSampleDir+ProjectPath.sampleID
                +"_control_Recal.bam -nct "+ProjectPath.thread)
        //MuTect2
        rec.execute("java -Xms20g -jar "+ProjectPath.GATK+" -T MuTect2 -R "+ProjectPath.hg38Reference+" -I:tumor "
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer_Recal.bam -I:normal "
                +ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_control_Recal.bam -o "
                +ProjectPath.sampleBaseDir+"somatic_mutations.vcf -nct "+ProjectPath.thread)
        //这个需要改最终的保存路径、其他一些参数

        //Vep注释之前需要先pass过滤

        //VEP
        rec.execute(ProjectPath.vep+" --cache --dir "+ProjectPath.vepReference+" --format vcf --input_file "+ProjectPath.sampleBaseDir
                +"somatic_mutations_pass_vepsource.vcf --output_file "+ProjectPath.sampleBaseDir+"somatic_mutations_pass_vep.vcf")

    }*/


}