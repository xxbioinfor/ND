package org.synergylab.neoantigenDiscovery.Reporter

import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.DeNovo.Trinity
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess.*
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess.CleanFile
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess.RawFile
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction.bindingPrediction
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction.NSSomaticMutation
import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.Prediction.generateHlaType
import org.synergylab.neoantigenDiscovery.utils.Filter
import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.getFileLines

//import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PreProcess.RawFile

fun main(args: Array<String>) {
    println("***** Start Discovering Neoantigens *****")

    println("===== Process Raw RNA-seq Files =====")
    val rawCancerSequence1 = ProjectPath.rawCancerSequence_1
    val rawCancerSequence2 = ProjectPath.rawCancerSequence_2
    val sampleID = ProjectPath.sampleID
    val cleanCancerSampleDir = ProjectPath.cleanCancerSampleDir
    val cancerSample = "cancer"
    val controlSample = "control"
    val cleanCancerFile = RawFile().trimmomaticCommand(rawCancerSequence1,rawCancerSequence2,cancerSample,sampleID,cleanCancerSampleDir)
    val cleanCancerSequence1 = cleanCancerFile.split("\t").get(0)
    val cleanCancerSequence2 = cleanCancerFile.split("\t").get(1)
    val fastqcCancerSequence1 = RawFile().fastqcCommand(cleanCancerSequence1,cleanCancerSampleDir)
    val fastqcCancerSequence2 = RawFile().fastqcCommand(cleanCancerSequence2,cleanCancerSampleDir)

    val rawControlSequence1 = ProjectPath.rawControlSequence_1
    val rawControlSequence2 = ProjectPath.rawControlSequence_2
    val cleanControlSampleDir = ProjectPath.cleanControlSampleDir
    val cleanControlFile = RawFile().trimmomaticCommand(rawControlSequence1,rawControlSequence2,controlSample,sampleID,cleanControlSampleDir)
    val cleanControlSequence1 = cleanControlFile.split("\t").get(0)
    val cleanControlSequence2 = cleanControlFile.split("\t").get(1)
    val fastqcControlSequence1 = RawFile().fastqcCommand(cleanControlSequence1,cleanControlSampleDir)
    val fastqcControlSequence2 = RawFile().fastqcCommand(cleanControlSequence2,cleanControlSampleDir)

    //generateCleanSequenceData()
    println(" Generate Clean RNA-seq Files! ")
    println(" Process Clean RNA-seq Files... ")
    val cancerSampleDir = ProjectPath.cancerSampleDir
    val cancerSamFile = CleanFile().generateSamFile(cleanCancerSequence1,cleanCancerSequence2,sampleID,cancerSampleDir)
    println(" Generate Cancer Sam File! ")
    val cancerBamFile = CleanFile().generateBamFile(cancerSamFile,sampleID,cancerSampleDir)
    println(" Generate Cancer Bam File! ")
    val controlSampleDir = ProjectPath.controlSampleDir
    val controlSamFile = CleanFile().generateSamFile(cleanControlSequence1,cleanControlSequence2,sampleID,controlSampleDir)
    println(" Generate Control Sam File! ")
    val controlBamFile = CleanFile().generateBamFile(controlSamFile,sampleID,controlSampleDir)
    println(" Generate Control Bam File! ")

    println("===== Process Clean RNA-seq Files =====")
    println(" Generate Somatic Mutation... ")
    val cancerPicardDir = ProjectPath.cancerSamplePicardDir
    val cancerReorderBamFile = NSSomaticMutation().picardCommand(cancerBamFile,cancerSample,sampleID,cancerPicardDir)
    val cancerGatkDir = ProjectPath.cancerSampleGatkDir
    val cancerRecalBamFile = NSSomaticMutation().gatkCommand(cancerReorderBamFile,sampleID,cancerGatkDir)
    val controlPicardDir = ProjectPath.controlSamplePicardDir
    val controlReorderBamFile = NSSomaticMutation().picardCommand(controlBamFile,controlSample,sampleID,controlPicardDir)
    val controlGatkDir = ProjectPath.controlSampleGatkDir
    val controlRecalBamFile = NSSomaticMutation().gatkCommand(controlReorderBamFile,sampleID,controlGatkDir)
    val mutationDir = ProjectPath.mutationDir
    val somaticMutation = NSSomaticMutation().mutect2Command(cancerRecalBamFile,controlRecalBamFile,mutationDir)
    val somaticMutationPass = Filter().filterGatkPass(somaticMutation,mutationDir)
    val mutationVepAnnotation = NSSomaticMutation().vepCommand(somaticMutationPass,mutationDir)
    println(" Generate Somatic Mutation! ")
    println(" Generate Trinity Transcript... ")
    val cancerTrinityFasta = Trinity().DeNovoTranscript(cleanCancerSequence1,cleanCancerSequence2,cancerSampleDir)
    val controlTrinityFasta = Trinity().DeNovoTranscript(cleanControlSequence1,cleanControlSequence2,controlSampleDir)
    println(" Generate Trinity Transcript! ")
    println(" Generate Trinity Transcript Expression... ")
    val cancerTrinityExpression = Trinity().DeNovoExpression(cancerTrinityFasta,cleanCancerSequence1,cleanCancerSequence2,cancerSampleDir)
    val controlTrinityExpression = Trinity().DeNovoExpression(controlTrinityFasta,cleanControlSequence1,cleanControlSequence2,controlSampleDir)
    println(" Generate Trinity Transcript Expression! ")
    /*
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
    */

    println("===== PostProcess =====")
    println(" Generate Mutated Transcript Files... ")
    val gmapDir = ProjectPath.gmapDir
    val gmapFile = gmapReference(cancerTrinityFasta,gmapDir)
    val gmapReference = MutatedTranscript().generateGmapReference(gmapFile,gmapDir)
    val sampleDir = ProjectPath.sampleBaseDir
    val mutatedTranscript = MutatedTranscript().generateMutatedTranscript(gmapReference,mutationVepAnnotation,sampleDir)
    println(" Generate Mutated Transcript Files! ")
    println(" Filter Info... ")
    val fpkmFilter = Filter().filterRsemFpkm(mutatedTranscript,cancerTrinityExpression,1.0, sampleDir)
    //val peptideCountsFilter = Filter().filterPeptideCounts(fpkmFilter,"",1.0,sampleDir)
    println(" Filter Info! ")
    println(" Generate HLA Type... ")
    val hlaFile = generateHlaType(cleanCancerSequence1,cleanCancerSequence2,sampleDir)
    val hlaType1 = getFileLines(hlaFile).get(1).split("\t").get(1)  //文件提取
    val hlaType2 = getFileLines(hlaFile).get(1).split("\t").get(3)
    var hlaType = ""
    if (hlaType1.equals(hlaType2)) {
        hlaType = hlaType1
    }
    else hlaType = hlaType1+","+hlaType2
    println(" Generate HLA Type: "+hlaType+"! ")
    for (i in 8..11) {
        println(" Generate "+i+"-mers Peptides Fasta File... ")
        val pepDir = ProjectPath.peptideDir
        val peptideFasta = PeptideFasta().generatePeptideFasta(fpkmFilter, i, pepDir)
        println(" Generate "+i+"-mers Peptides Fasta File! ")
        println(" Predict "+i+"-mers Peptide／MHC Binding Affinities...")
        val affinityPrediction = bindingPrediction(peptideFasta,hlaType, sampleID, pepDir)
        println(" Predict "+i+"-mers Peptide／MHC Binding Affinities! ")
        val affinityFilter = Filter().filterBindingAffinity(affinityPrediction,hlaType,pepDir)
    }
    val neoantigen = ""

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
