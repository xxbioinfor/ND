package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.DeNovo

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//trinity从头拼接的转录本表达量信息

//class DeNovoTranscriptExpression(){
    //parameter

    //source data
    //val denovoTranscript = 0 //来自cleanfileprocess处理得到的bamfile

    //class properties
    //val denovoTranscriptExpression ＝ generateDeNovoTranscriptExpression()

    //function
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

//}