package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//预测peptide与MHC结合亲和力

//class pMHCBinding(){
    //parameters

    //source data
    //val peptideFasta = 0 //读取generatepeptidefasta文件
    //val hlaType = 0

    //class properties
    //val bindingResult = bindingPrediction()

    //function
    fun bindingPrediction(peptidesFastaFile: String,hlaType: String,sampleID: String,outDir: String): String {
        //使用NETMHC 4.0进行预测
        //将得到的结果文件进行整理，并按bindingresult的格式输出
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP, ProjectPath.linuxName, ProjectPath.linuxPwd)
        val netMHC = ProjectPath.netMHC
        rec.execute("${netMHC} -a ${hlaType} -xls -f ${peptidesFastaFile} -xlsfile ${outDir}${sampleID}_binding.xls")

        return "${outDir}${sampleID}_binding.xls"

    }

//}