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
    fun bindingPrediction(){
        //使用NETMHC 4.0进行预测
        //将得到的结果文件进行整理，并按bindingresult的格式输出
        val hlaType = ""
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP, ProjectPath.linuxName, ProjectPath.linuxPwd)
        rec.execute(ProjectPath.netMHC+" -a "+hlaType+" -xls -f "+ProjectPath.cancerSampleDir+ProjectPath.sampleID
                +"_peptides.fsa -xlsfile "+ProjectPath.cancerSampleDir+ProjectPath.sampleID+"_binding.xls")

    }

//}