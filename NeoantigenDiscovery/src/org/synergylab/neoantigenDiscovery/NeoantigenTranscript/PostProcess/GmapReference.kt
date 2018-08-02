package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//class GmapReference(){

    //function
    fun gmapReference(){
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP, ProjectPath.linuxName, ProjectPath.linuxPwd)
        rec.execute("gmap -D "+ProjectPath.gmap+" -d Homo_sapiens_assembly38 -Z -t "+ProjectPath.thread
                +ProjectPath.cancerSampleDir+"trinity/Trinity.fasta > "+ProjectPath.cancerSampleDir+"gmap_out")
    }
//}