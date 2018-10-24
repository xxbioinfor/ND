package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//由bam格式文件,使用seq2HLA预测HLA分型

//class HLAType(){
    //parameter

    //source data
    //val bamFile = 0 //来自cleanfileprocess处理得到的bamfile

    //class properties
    //val hlaType ＝ generateHlaType()

    //function
    fun generateHlaType(cleanFile1: String,cleanFile2: String,outDir: String): String {
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP, ProjectPath.linuxName, ProjectPath.linuxPwd)
        //rec.execute("mkdir "+ProjectPath.cancerSampleDir+"hla")
        rec.execute("python2.7 "+ProjectPath.seq2HLA+" -1 "+cleanFile1+" -2 "+cleanFile2+" -r "
                +outDir+"hla/hla -p "+ProjectPath.thread)
        //rec.execute("python2.7 /localdisk/software/seq2HLA2.2/seq2HLA.py -1 /localdisk/jobs/lung/180313002ML/cancer/cleandata/180313002ML_cancer_fq1_clean.fq "
        //        +"-2 /localdisk/jobs/lung/180313002ML/cancer/cleandata/180313002ML_cancer_fq2_clean.fq -r /localdisk/jobs/lung/180313002ML/cancer/hla/180313002ML -p 10")
        //注意这个要用python2版本运行
        return outDir+"hla-ClassI.HLAgenotype4digits"
    }

//}