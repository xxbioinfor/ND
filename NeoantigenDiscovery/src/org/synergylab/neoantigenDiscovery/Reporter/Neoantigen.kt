package org.synergylab.neoantigenDiscovery.Reporter

import NeoantigenDiscovery.NeoantigenResult
import org.synergylab.neoantigenDiscovery.NeoantigenMS.generateMzxmlFile
import org.synergylab.neoantigenDiscovery.NeoantigenMS.generatePeptideCounts
import org.synergylab.neoantigenDiscovery.NeoantigenMS.generatePeptideTranscriptFile
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.*
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.DeNovo.generateDeNovoTranscript
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.DeNovo.generateDeNovoTranscriptExpression
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess.generateGmapReference
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess.generateMutatedTranscript
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess.generatePeptideFasta
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess.generateTranscriptTranslation
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess.generateBamFile
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess.generateCleanSequenceData
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess.generateSamFile
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction.bindingPrediction
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction.generateHlaType
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction.generateSomaticMutation
import org.synergylab.neoantigenDiscovery.utils.filterFPKM
import org.synergylab.neoantigenDiscovery.utils.filterPeptideCounts

//import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess.RawFile

fun main(args: Array<String>) {
    println("***** Start Discovering Neoantigens *****")

    println("===== Process Raw RNA-seq Files =====")
    generateCleanSequenceData()
    println(" Generate Clean RNA-seq Files! ")
    println(" Process Clean RNA-seq Files... ")
    generateSamFile()
    println(" Generate Sam Files! ")
    generateBamFile()
    println(" Generate Bam Files! ")

    println("===== Process Clean RNA-seq Files =====")
    println(" Generate Somatic Mutation... ")
    generateSomaticMutation()
    println(" Generate Somatic Mutation! ")
    println(" Generate Trinity Transcript... ")
    generateDeNovoTranscript()
    println(" Generate Trinity Transcript! ")
    println(" Generate Trinity Transcript Expression... ")
    generateDeNovoTranscriptExpression()
    println(" Generate Trinity Transcript Expression! ")
    println(" Generate Protein Database... ")
    generateTranscriptTranslation()
    println(" Generate Protein Database! ")

    println("===== Process MS Files =====")
    println(" Generate MzXML Files... ")
    generateMzxmlFile()
    println(" Generate MzXML Files! ")
    println(" Generate Peptide Counts Files... ")
    generatePeptideTranscriptFile()
    generatePeptideCounts()
    println(" Generate Peptide Counts Files! ")

    println("===== PostProcess =====")
    println(" Generate Mutated Transcript Files... ")
    generateGmapReference()
    generateMutatedTranscript()
    println(" Generate Mutated Transcript Files! ")
    println(" Filter Info... ")
    filterFPKM()
    filterPeptideCounts()
    println(" Filter Info! ")
    println(" Generate HLA Type... ")
    generateHlaType()
    println(" Generate HLA Type! ")
    println(" Generate Peptides Fasta Files... ")
    generatePeptideFasta()
    println(" Generate Peptides Fasta Files! ")
    println(" Predict Binding Affinity... ")
    bindingPrediction()
    println(" Predict Binding Affinity! ")
    //filter()




}

/*
class Neoantigen(sampleName: String){
    //parameters
    val sample = sampleName

    //source data
    val pMHCFile = 0 //读取pMHCBinding相应文件
    RawFile()
    CleanFile()
    DeNovoAssembly()


    //class properties
    val neoantigen = generateNeoantigenResult()

    //function
    private fun generateNeoantigenResult(): NeoantigenResult{
        //提取相应信息
        //以文件形式输出
    }
}
*/
