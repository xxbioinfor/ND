package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess

import org.synergylab.neoantigenDiscovery.utils.appendFile
import org.synergylab.neoantigenDiscovery.utils.getFileLines
import org.synergylab.neoantigenDiscovery.utils.getSubUtilSimple


//生成(2x－1)AAs肽段
//class PeptideFasta(){
    //parameters

    //source data
    //val mutation = 0 //筛选出的突变信息文件
    //val mutationSite = 0 //突变信息中的突变位点
    //val mutationType = 0 //突变信息中的突变类型

    //class properties
    //val peptideFasta = generatePeptideFasta()

    //function
    fun generatePeptideFasta(): String {
        val mutatedTranscriptFileName = "/Users/toby/Desktop/neoantigenData/MutatedTranscript_pepCounts1.0.txt"
        val mutatedTranscriptFile = getFileLines(mutatedTranscriptFileName)
        val mutatedTranscriptSize = mutatedTranscriptFile.size

        val peptideReferenceFileName = "/Users/toby/Desktop/neoantigenData/Homo_sapiens.GRCh38.pep.all.fa"
        val peptideReferenceFile = getFileLines(peptideReferenceFileName)
        val peptideReferenceSize = peptideReferenceFile.size

        val peptideReferenceMap = HashMap<String,String>()
        var transcriptID = ""
        for (i in 0 until peptideReferenceSize){
            val peptideReferenceLine = peptideReferenceFile.get(i)
            if (peptideReferenceLine.startsWith(">")){
                transcriptID = getSubUtilSimple(peptideReferenceLine,"transcript:(.*)\\.")
            }
            else peptideReferenceMap.replace(transcriptID,peptideReferenceMap[transcriptID]+peptideReferenceLine)
        }

        val peptideFastaLength = 21
        var mutatedPeptide = ""
        var wildPosition = 0
        var wildPeptide = ""
        val peptideFastaFileName = "/Users/toby/Desktop/neoantigenData/peptideFasta.fasta"

        for (i in 0 until mutatedTranscriptSize) {
            val mutatedTranscriptLine = mutatedTranscriptFile.get(i)
            val mtTrinityID = mutatedTranscriptLine.split("\t").get(0)
            //val chr = mutatedTranscriptLine.split("\t").get(1)
            //val alt = mutatedTranscriptLine.split("\t").get(2)
            //val mtStrand = mutatedTranscriptLine.split("\t").get(3)
            //val trinity_start = mutatedTranscriptLine.split("\t").get(4).toInt()
            //val trinity_end = mutatedTranscriptLine.split("\t").get(5).toInt()
            //val ref_start = mutatedTranscriptLine.split("\t").get(6).toInt()
            //val ref_end = mutatedTranscriptLine.split("\t").get(7).toInt()
            //val ref_position_start = mutatedTranscriptLine.split("\t").get(8).toInt()
            //val ref_position_end = mutatedTranscriptLine.split("\t").get(9).toInt()
            val transcript = mutatedTranscriptLine.split("\t").get(10)
            val mutationType = mutatedTranscriptLine.split("\t").get(11)
            val proteinPosition = mutatedTranscriptLine.split("\t").get(12).toInt()
            val aminoAcid = mutatedTranscriptLine.split("\t").get(13).split("/").get(1)

            val wildPeptideResult = generateWildSequence(proteinPosition, peptideReferenceMap[transcript].toString(),peptideFastaLength)
            wildPeptide = wildPeptideResult.split("\t").get(0)
            wildPosition = wildPeptideResult.split("\t").get(1).toInt()
            mutatedPeptide = generateMutatedSequence(wildPosition,wildPeptide,aminoAcid,mutationType)

            if (!wildPeptide.equals("") && !mutatedPeptide.equals("")){
                val wildHeader = ">wild "+mtTrinityID+" "+transcript+"\n"
                val mutatedHeader = ">mutated "+mtTrinityID+" "+transcript+"\n"
                appendFile(wildHeader+wildPeptide+"\n"+mutatedHeader+mutatedPeptide+"\n",peptideFastaFileName)
            }
        }

        return peptideFastaFileName
    }


    private fun generateDistanceFromStart(position: Int, sequence: String): Int {
        return position
    }
    private fun generateDistanceFromEnd(position: Int, sequence: String): Int{
        return sequence.length-1-position
    }
    private fun generateFlankingLength(peptideFastaLength: Int): Int {
        if (peptideFastaLength % 2 == 0)
            return (peptideFastaLength-2)/2
        else
            return (peptideFastaLength-1)/2
    }

    private fun generateWildSequence(position: Int, sequence: String, peptideFastaLength: Int): String {
        //position为蛋白突变位点
        //sequence为蛋白序列
        //peptideFastaLength为肽段全长
        var wildPeptide = ""
        var wildPosition = 0
        var wildResult = ""

        val distanceFromStart = generateDistanceFromStart(position,sequence)
        val distanceFromEnd = generateDistanceFromEnd(position,sequence)
        val flankingLength = generateFlankingLength(peptideFastaLength)
        val sequenceLength = sequence.length

        if (distanceFromStart < flankingLength){
            wildPeptide = sequence.substring(0..peptideFastaLength)
            wildPosition = position //突变位于肽段的位置
            wildResult = wildPeptide+"\t"+wildPosition
        }
        else if (distanceFromEnd < flankingLength){
            val startPosition = sequenceLength - peptideFastaLength
            wildPeptide = sequence.substring(startPosition..sequenceLength)
            wildPosition = peptideFastaLength - distanceFromEnd - 1
            wildResult = wildPeptide+"\t"+wildPosition
        }
        else if (distanceFromStart >= flankingLength && distanceFromEnd >= flankingLength){
            val startPosition = position - flankingLength
            val endPosition = startPosition + peptideFastaLength
            wildPeptide = sequence.substring(startPosition,endPosition)
            wildPosition = flankingLength
            wildResult = wildPeptide+"\t"+wildPosition
        }
        return wildResult

    }

    private fun generateMutatedSequence(wildPosition: Int,wildPeptide: String,mutationBases: String,mutationType: String): String {
        //position是指wildPosition
        var mutatedPeptide = ""

        val mutatedBase = mutationBases.split("/").get(1)
        val BaseNum = mutatedBase.length
        val sequenceLength = mutatedPeptide.length
        val mutationType = mutationType
        val position = wildPosition
        val wildPeptide = wildPeptide

        if (mutationType.equals("missense_variant")){
            mutatedPeptide = wildPeptide.substring(0..position)+mutatedBase+wildPeptide.substring(position+1..sequenceLength)
        }
        if (mutationType.equals("inframe_deletion")){
            mutatedPeptide = wildPeptide.substring(0..position)+wildPeptide.substring(position+1..sequenceLength+BaseNum)
        }
        if (mutationType.equals("inframe_insertion")){
            mutatedPeptide = wildPeptide.substring(0..position)+mutatedBase+wildPeptide.substring(position+1..sequenceLength-BaseNum)
        }
        return mutatedPeptide
    }

    /*
    fun generatePeptideFasta(){
        //由突变位点、突变类型信息
        //生成peptidefasta
        val mutatedTranscriptFileName = "/Users/toby/Desktop/neoantigenData/MutatedTranscript_pepCounts1.0.txt"
        val mutatedTranscriptFile = getFileLines(mutatedTranscriptFileName)
        val mutatedTranscriptSize = mutatedTranscriptFile.size

        //trinity_cancer
        val trinityCancerFileName = "/Users/toby/Desktop/neoantigenData/trinity_cancer"
        val trinityCancerFile = getFileLines(trinityCancerFileName)
        val trinityCancerSize = trinityCancerFile.size
        /*
        val trinityCancerMap = HashMap<String,String>()
        for (i in 0 until trinityCancerSize){
            val trinityCancerLine = trinityCancerFile.get(i)
            val index = trinityCancerLine.split("\t").get(0)
            trinityCancerMap[index] = trinityCancerLine
        }
        */

        //trinity_protein_sequence
        val proteinSequenceFileName = "/Users/toby/Desktop/neoantigenData/trinity_cancer_fasta.tab"
        val proteinSequenceFile = getFileLines(proteinSequenceFileName)
        val proteinSequenceSize = proteinSequenceFile.size
        val proteinSequenceMap = HashMap<String,String>()
        var proteinIndex = ""
        for (i in 0 until proteinSequenceSize){
            val proteinSequenceLine = proteinSequenceFile.get(i)
            if (proteinSequenceLine.startsWith(">")){
                proteinIndex = proteinSequenceLine.substring(1).split("|").get(0)
            }else{
                proteinSequenceMap[proteinIndex] = proteinSequenceLine
            }
        }

        var mutatedPeptideResult = ""
        var mutationPositionResult = 0
        var wildPeptideResult = ""

        for (i in 0 until mutatedTranscriptSize){
            val mutatedTranscriptLine = mutatedTranscriptFile.get(i)
            val mtTrinityID = mutatedTranscriptLine.split("\t").get(0)
            val chr = mutatedTranscriptLine.split("\t").get(1)
            val alt = mutatedTranscriptLine.split("\t").get(2)
            val mtStrand = mutatedTranscriptLine.split("\t").get(3)
            val trinity_start = mutatedTranscriptLine.split("\t").get(4).toInt()
            val trinity_end = mutatedTranscriptLine.split("\t").get(5).toInt()
            val ref_start = mutatedTranscriptLine.split("\t").get(6).toInt()
            val ref_end = mutatedTranscriptLine.split("\t").get(7).toInt()
            val ref_position_start = mutatedTranscriptLine.split("\t").get(8).toInt()
            val ref_position_end = mutatedTranscriptLine.split("\t").get(9).toInt()
            var trinity_position_start = 0
            var trinity_position_end = "-"

            if (mtStrand.equals("+")){
                trinity_position_start = trinity_start + (ref_position_start - ref_start)
                if (!ref_position_end.equals("-")){
                    trinity_position_end = (trinity_start + (ref_position_start - ref_start)).toString()
                }
            }else if (mtStrand.equals("-")){
                trinity_position_start = trinity_start + (ref_position_start - ref_end)
                if (!ref_position_end.equals("-")){
                    trinity_position_end = (trinity_start + (ref_position_start - ref_end)).toString()
                }
            }

            for (j in 0 until trinityCancerSize){
                val trinityCancerLine = trinityCancerFile.get(i)
                val tcTrinityID = trinityCancerLine.split("\t").get(1)
                if (tcTrinityID.equals(mtTrinityID)){
                    val index = trinityCancerLine.split("\t").get(0)
                    val tcStrand = trinityCancerLine.split("\t").get(2)
                    val frame = trinityCancerLine.split("\t").get(3)
                    val tc_start = trinityCancerLine.split("\t").get(4).toInt()
                    val tc_end = trinityCancerLine.split("\t").get(5).toInt()
                    var mutatedSequenceResult = ""

                    if (trinity_position_end.equals("-") && tc_start <= trinity_position_start && trinity_position_start <= tc_end){
                        mutatedSequenceResult = generateWildSequence(trinity_position_start,proteinSequenceMap[index].toString(),21)
                        mutatedPeptideResult = mutatedSequenceResult.split("\t").get(0)
                        mutationPositionResult = mutatedSequenceResult.split("\t").get(1).toInt()
                    }
                    if (!trinity_position_end.equals("-") && tc_start <= trinity_position_start && trinity_position_end.toInt() <= tc_end){

                    }
                }




            }

        }

    }*/
    /*
    private fun generateMutatedSequence(position: Int, sequence: String, peptideFastaLength: Int): String {
        var mutatedPeptide = ""
        var mutationPosition = 0
        //val resultMap = HashMap<String,Int>()
        var mutatedResult = ""

        val distanceFromStart = generateDistanceFromStart(position,sequence)
        val distanceFromEnd = generateDistanceFromEnd(position,sequence)
        val flankingLength = generateFlankingLength(peptideFastaLength)
        val sequenceLength = sequence.length
        if (distanceFromStart < flankingLength){
            mutatedPeptide = sequence.substring(0..peptideFastaLength)
            mutationPosition = position //突变位于突变肽段的位置
            mutatedResult = mutatedPeptide+"\t"+mutationPosition
            //if (!resultMap.containsKey(mutatedPeptide)) {
            //    resultMap.put(mutatedPeptide, mutationPosition)
            //}
        }
        else if (distanceFromEnd < flankingLength){
            val startPosition = sequenceLength - peptideFastaLength
            mutatedPeptide = sequence.substring(startPosition..sequenceLength)
            mutationPosition = peptideFastaLength - distanceFromEnd - 1
            mutatedResult = mutatedPeptide+"\t"+mutationPosition
            //if (!resultMap.containsKey(mutatedPeptide)) {
            //    resultMap.put(mutatedPeptide, mutationPosition)
            //}
        }
        else if (distanceFromStart >= flankingLength && distanceFromEnd >= flankingLength){
            val startPosition = position - flankingLength
            val endPosition = startPosition + peptideFastaLength
            mutatedPeptide = sequence.substring(startPosition,endPosition)
            mutationPosition = flankingLength
            mutatedResult = mutatedPeptide+"\t"+mutationPosition
            //if (!resultMap.containsKey(mutatedPeptide)) {
            //    resultMap.put(mutatedPeptide, mutationPosition)
            //}
        }
        //return resultMap //返回两个？ mutationPosition
        return mutatedResult
    }
    private fun generateWildPeptide(mutatedPeptide: String,mutationPosition: Int,mutationBases: String,mutationType: String): String {
        //直接替换 突变前序列+突变+突变后序列
        val wildBase = mutationBases.split("/").get(1) //先暂定前面的为野生
        val sequenceLength = mutatedPeptide.length
        val mutationType = mutationType
        var wildPeptide = ""
        if ("missense_variant".toRegex().containsMatchIn(mutationType) || "inframe_deletion".toRegex().containsMatchIn(mutationType) || "inframe_insertion".toRegex().containsMatchIn(mutationType)){
            wildPeptide = mutatedPeptide.substring(0..mutationPosition)+wildBase+mutatedPeptide.substring(mutationPosition+2..sequenceLength)
        }
        //if ("frameshift_variant".toRegex().containsMatchIn(mutationType)){
        //    val wildPeptide = mutatedPeptide.substring(0..mutationPosition)
        //}
        return wildPeptide

    }
    */
//}