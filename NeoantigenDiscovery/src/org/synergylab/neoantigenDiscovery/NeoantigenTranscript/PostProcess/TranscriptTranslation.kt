package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//对从头拼接数据，进行六框翻译

//class TranscriptTranslation(){
    //parameter

    //source data
    //val denovoTranscript = 0

    //class properties
    //val transcriptTranslation ＝ generateTranscriptTranslation()

    //function
    fun generateTranscriptTranslation(){
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)
        rec.execute("R")
        rec.execute("library('PGA')")
        rec.execute("transcript_seq_file <- \""+ProjectPath.cancerSampleDir+"trinity/Trinity.fasta\"")
        rec.execute("outfile_name <- \""+ProjectPath.cancerSampleDir+"trinity/trinity_cancer\"")
        //结果文件为 trinity_cancer_txFinder.fasta
        rec.execute("outmtab <- \""+ProjectPath.cancerSampleDir+"trinity/trinity_cancer_fasta.tab\"")
        rec.execute("createProDB4DenovoRNASeq(infa=transcript_seq_file,outfile_name,outmtab, bool_use_3frame=FALSE, bool_get_longest=FALSE)")
        rec.execute("q()")
        rec.execute("n")
    }

//}