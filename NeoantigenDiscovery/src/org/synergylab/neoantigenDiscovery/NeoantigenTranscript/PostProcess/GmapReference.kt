package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//class GmapReference(){

    //function
    fun gmapReference(trinityFile: String,outDir: String): String {
        val rec = RemoteExecuteCommand(ProjectPath.linuxIP, ProjectPath.linuxName, ProjectPath.linuxPwd)
        val gmap = ProjectPath.gmap
        val thread = ProjectPath.thread
        rec.execute("gmap -D ${gmap} -d Homo_sapiens_assembly38 -A -t ${thread} ${trinityFile} > ${outDir}gmap_out")
        return "${outDir}gmap_out"
        //-A -Z
    }
//}