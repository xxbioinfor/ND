package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.DeNovo

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//对clean.fq数据
//使用triniy进行从头组装

class Trinity() {
    val rec = RemoteExecuteCommand(ProjectPath.linuxIP, ProjectPath.linuxName, ProjectPath.linuxPwd)
    val trinity = ProjectPath.trinity
    val thread = ProjectPath.thread
    fun DeNovoTranscript(cleanFile1: String, cleanFile2: String, outDir: String): String {
        rec.execute("${trinity}Trinity --seqType fq --max_memory 20G --left ${cleanFile1} --right ${cleanFile2}" +
                " --CPU ${thread} --output ${outDir}trinity")
        return "${outDir}trinity/Trinity.fasta"
    }

    fun DeNovoExpression(trinityFile: String,cleanFile1: String,cleanFile2: String,outDir: String): String {
        rec.execute("${trinity}util/align_and_estimate_abundance.pl --transcript ${trinityFile} --seqType fq " +
                "--left ${cleanFile1} --right ${cleanFile2} --est_method RSEM --aln_method bowtie2 --prep_reference " +
                "--trinity_mode --output_dir ${outDir}trinity")
        return "${outDir}trinity/RSEM.isoforms.results"
    }

/*
    //function
    fun generateDeNovoTranscript(){
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)
        //cancer
        //rec.execute("/localdisk/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --max_memory 20G " +
        //        "--left /localdisk/jobs/lung/180313002ML/cancer/cleandata/180313002ML_cancer_fq1_clean.fq " +
        //        "--right /localdisk/jobs/lung/180313002ML/cancer/cleandata/180313002ML_cancer_fq2_clean.fq --CPU 20 " +
        //        "--output /localdisk/jobs/lung/180313002ML/cancer/trinity_out_dir")
        rec.execute(ProjectPath.trinity+"Trinity --seqType fq --max_memory 20G --left "+ProjectPath.cleanCancerSampleDir
                +ProjectPath.sampleID+"_cancer_fq1_clean.fq --right "+ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID
                +"_cancer_fq2_clean.fq --CPU "+ProjectPath.thread+" --output "+ProjectPath.cancerSampleDir+"trinity_out_dir")
        //control
        //rec.execute("/localdisk/software/trinityrnaseq-Trinity-v2.5.1/Trinity --seqType fq --max_memory 20G " +
        //        "--left /localdisk/jobs/lung/180313002ML/control/cleandata/180313002ML_control_fq1_clean.fq " +
        //        "--right /localdisk/jobs/lung/180313002ML/control/cleandata/180313002ML_control_fq2_clean.fq --CPU 20 " +
        //        "--output /localdisk/jobs/lung/180313002ML/control/trinity_out_dir")
        rec.execute(ProjectPath.trinity+"Trinity --seqType fq --max_memory 20G --left "+ProjectPath.cleanControlSampleDir
                +ProjectPath.sampleID +"_control_fq1_clean.fq --right "+ProjectPath.cleanControlSampleDir+ProjectPath.sampleID
                +"_control_fq2_clean.fq --CPU "+ProjectPath.thread+" --output "+ProjectPath.controlSampleDir+"trinity_out_dir")
    }


    fun generateDeNovoTranscriptExpression(){
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)

        //cancer
        rec.execute(ProjectPath.trinity+"util/align_and_estimate_abundance.pl --transcript "+ProjectPath.cancerSampleDir
                +"trinity/Trinity.fasta --seqType fq --left "+ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID
                +"_cancer_fq1_clean.fq --right "+ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID
                +"_cancer_fq2_clean.fq --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode --output_dir "
                +ProjectPath.cancerSampleDir+"trinity")
        //rec.execute("/localdisk/software/trinityrnaseq-Trinity-v2.5.1/util/align_and_estimate_abundance.pl --transcript /localdisk/jobs/lung/180313002ML/cancer/trinity/Trinity.fasta --seqType fq --left /localdisk/jobs/lung/180313002ML/cancer/cleandata/180313002ML_cancer_fq1_clean.fq --right /localdisk/jobs/lung/180313002ML/cancer/cleandata/180313002ML_cancer_fq2_clean.fq --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode --output_dir /localdisk/jobs/lung/180313002ML/cancer/trinity")

        //control
        rec.execute(ProjectPath.trinity+"util/align_and_estimate_abundance.pl --transcript "+ProjectPath.controlSampleDir
                +"trinity/Trinity.fasta --seqType fq --left "+ProjectPath.cleanControlSampleDir+ProjectPath.sampleID
                +"_control_fq1_clean.fq --right "+ProjectPath.cleanControlSampleDir+ProjectPath.sampleID
                +"_control_fq2_clean.fq --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode --output_dir "
                +ProjectPath.controlSampleDir+"trinity")
    }
*/
}

