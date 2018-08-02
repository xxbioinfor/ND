package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//将clean.fq数据处理成bam格式

//class CleanFile(){

    //val samFile = generateSamFile()
    //val bamFile = generateBamFile()

    //function
    fun generateSamFile(){
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)
        //hisat2
        //cancer
        rec.execute(ProjectPath.hisat2+" -p 20 -x "+ProjectPath.hisat2Index+" -1 "+ProjectPath.cleanCancerSampleDir
                +ProjectPath.sampleID+"_cancer_fq1_clean.fq -2 "+ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID
                +"_cancer_fq2_clean.fq -S "+ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer.sam")
        //control
        rec.execute(ProjectPath.hisat2+" -p 20 -x "+ProjectPath.hisat2Index+" -1 "+ProjectPath.cleanControlSampleDir
                +ProjectPath.sampleID+"_control_fq1_clean.fq -2 "+ProjectPath.cleanControlSampleDir+ProjectPath.sampleID
                +"_control_fq2_clean.fq -S "+ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control.sam")

    }

    fun generateBamFile(){
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)
        //samtools
        //cancer
        rec.execute(ProjectPath.samtools+" view -S -b "+ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer.sam | "
                +ProjectPath.samtools+" sort | "+ProjectPath.samtools+" view -h | perl -ne 'if(/HI:i:(\\d+)/) { \$m=\$m1-1; \$_ =~ s/HI:i:(\\d+)/HI:i:\$m/} print \$_;' | "
                +ProjectPath.samtools+" view -bS - > "+ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_cancer.bam")
        //control
        rec.execute(ProjectPath.samtools+" view -S -b "+ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control.sam | "
                +ProjectPath.samtools+" sort | "+ProjectPath.samtools+" view -h | perl -ne 'if(/HI:i:(\\d+)/) { \$m=\$m1-1; \$_ =~ s/HI:i:(\\d+)/HI:i:\$m/} print \$_;' | "
                +ProjectPath.samtools+" view -bS - > "+ProjectPath.controlSampleDir+ProjectPath.sampleID+"_control.bam")
    }


//}