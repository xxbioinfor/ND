package org.synergylab.neoantigenDiscovery.test

import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.DeNovo.generateDeNovoTranscript
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.DeNovo.generateDeNovoTranscriptExpression
import org.synergylab.neoantigenDiscovery.utils.ProjectPath

fun main(args: Array<String>) {
    //println(" Generate DeNovo Transcript... ")
    //try {
    //    generateDeNovoTranscript()
    //}catch (e: Exception){
    //    println(e.message)
    //}
    //println(" Generate DeNovo Transcript! ")
    println(" Generate DeNovo Transcript Expression... ")
    generateDeNovoTranscriptExpression()
    println(" Generate DeNovo Transcript Expression! ")
    //println(ProjectPath.trinity)
    //println(ProjectPath.cleanCancerSampleDir +ProjectPath.sampleID +"_cancer_fq1_clean.fq")
    //println(ProjectPath.cleanCancerSampleDir+ProjectPath.sampleID +"_cancer_fq2_clean.fq")
    //println(ProjectPath.thread)
    //println(ProjectPath.cancerSampleDir+"trinity_out_dir")
}