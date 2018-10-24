package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//预处理raw数据
//将raw数据进行质控、去adapter、去掉前15bp等操作，生成clean.fq

class RawFile(){
    val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)

    fun trimmomaticCommand(rawFile1: String,rawFile2: String,sampleName: String,sampleID: String,cleanSampleDir: String): String {
        val rawSequenceFile1 = rawFile1
        val rawSequenceFile2 = rawFile2
        val sample = sampleName
        val sampleID = sampleID
        val cleanSampleDir = cleanSampleDir

        rec.execute("java -jar "+ProjectPath.trimmomatic+" PE -threads "+ProjectPath.thread+" -phred33 "+rawSequenceFile1
                +" "+rawSequenceFile2+" "+cleanSampleDir+sampleID +"_"+sample+"_fq1_clean.fq "
                +cleanSampleDir+sampleID+"_"+sample+"_fq1_unpaired.fq "+cleanSampleDir+sampleID+"_"+sample
                +"_fq2_clean.fq "+cleanSampleDir +sampleID+"_"+sample+"_fq2_unpaired.fq ILLUMINACLIP:"
                +ProjectPath.adapter +":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 HEADCROP:15")

        val cleanFile = cleanSampleDir+sampleID +"_"+sample+"_fq1_clean.fq\t"+cleanSampleDir+sampleID +"_"+sample+"_fq2_clean.fq"
        return cleanFile
    }

    fun fastqcCommand(cleanFile: String,sampleDir: String): String {
        val cleanSequenceFile = cleanFile
        val cleanSampleDir = sampleDir
        rec.execute(ProjectPath.fastqc+" --extract -f fastq -t "+ProjectPath.thread+" "+cleanSequenceFile+" -o "+ProjectPath.cleanCancerSampleDir)
        val fileName = cleanSequenceFile.split(".").get(0)
        val fastqcFile = fileName+"_fastqc/fastqc_data.txt"  //文件可以再修改
        return fastqcFile

    }
    /*
    //function
    fun generateCleanSequenceData(){
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)

        //trimmomatic:去掉前15bp,去adapter,生成clean.fq
        //rec.execute("mkdir "+ProjectPath.sampleBaseDir+"cleandata") //创建cleandata文件夹
        //cancer
        rec.execute("java -jar "+ProjectPath.trimmomatic+" PE -threads "+ProjectPath.thread+" -phred33 "+ProjectPath.rawCancerSequence_1
                +" "+ProjectPath.rawCancerSequence_2+" "+ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID
                +"_cancer_fq1_clean.fq "+ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID+"_cancer_fq1_unpaired.fq "
                +ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID+"_cancer_fq2_clean.fq "+ProjectPath.cleanCancerSampleDir
                +ProjectPath.sampleID+"_cancer_fq2_unpaired.fq ILLUMINACLIP:"+ProjectPath.adapter
                +":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 HEADCROP:15")
        //control
        rec.execute("java -jar "+ProjectPath.trimmomatic+" PE -threads "+ProjectPath.thread+" -phred33 "+ProjectPath.rawControlSequence_1
                +" "+ProjectPath.rawControlSequence_2+" "+ProjectPath.cleanControlSampleDir+ProjectPath.sampleID
                +"_control_fq1_clean.fq "+ProjectPath.cleanControlSampleDir+ProjectPath.sampleID+"_control_fq1_unpaired.fq "
                +ProjectPath.cleanControlSampleDir+ProjectPath.sampleID+"_control_fq2_clean.fq "+ProjectPath.cleanControlSampleDir
                +ProjectPath.sampleID+"_control_fq2_unpaired.fq ILLUMINACLIP:"+ProjectPath.adapter
                +":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 HEADCROP:15")

        //fastqc:质控
        //cancer
        rec.execute(ProjectPath.fastqc+" --extract -f fastq -t "+ProjectPath.thread+" "+ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID
                +"_cancer_fq1_clean.fq " +"-o "+ProjectPath.cleanCancerSampleDir)
        rec.execute(ProjectPath.fastqc+" --extract -f fastq -t "+ProjectPath.thread+" "+ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID
                +"_cancer_fq2_clean.fq " +"-o "+ProjectPath.cleanCancerSampleDir)
        //control
        rec.execute(ProjectPath.fastqc+" --extract -f fastq -t "+ProjectPath.thread+" "+ProjectPath.cleanControlSampleDir+ProjectPath.sampleID
                +"_control_fq1_clean.fq " +"-o "+ProjectPath.cleanControlSampleDir)
        rec.execute(ProjectPath.fastqc+" --extract -f fastq -t "+ProjectPath.thread+" "+ProjectPath.cleanControlSampleDir+ProjectPath.sampleID
                +"_control_fq2_clean.fq " +"-o "+ProjectPath.cleanControlSampleDir)

    }
    */
}