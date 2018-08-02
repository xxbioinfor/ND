package org.synergylab.neoantigenDiscovery.NeoantigenMS

import org.synergylab.neoantigenDiscovery.utils.*
import java.io.File

//提取肽段与转录本信息
//得到肽段支持较多的转录本

//class PeptideCounts(){
    //parameters

    //source data
    //val pepXMLFile = 0 //读取pepxml文件

    //class properties
    //internal val peptideTranscriptFile = generatePeptideTranscriptFile()
    //internal val peptideCounts = generatePeptideCounts()

    //function
    fun generatePeptideTranscriptFile(): String {
        //val filename = ProjectPath.msDir+"*.tandem.pep.xml"
        val filename = "/Users/toby/Desktop/neoantigenData/20141110_QEp1_LC7_PhGe_4_53_M5_1.2018_07_10_12_43_54.t.tandem.pep.xml"
        var file = getFileLines(filename)
        var size = file.size
        //val peptideTranscriptFilePath = ProjectPath.msDir+"peptideTranscript.txt"
        val peptideTranscriptFilePath = "/Users/toby/Desktop/neoantigenData/peptideTranscript_2.txt"

        val re_spec = "<spectrum_query spectrum=\"(.*?)\"" 
        val re_search = "<search_hit"
        val re_pep = "peptide=\"(.*?)\""
        val re_pro = "protein=\"(.*?)\""
        val re_alter = "<alternative_protein protein=\"(.*?)\""

        var spectrum = ""
        var peptide = ""
        var protein = ""
        var transcript = ""
        var outLine = ""

        for(i in 0..size-1){
            var line = file.get(i)

            if (re_spec.toRegex().containsMatchIn(line)) {
                spectrum = getSubUtilSimple(line, re_spec)
            }
            if (re_search.toRegex().containsMatchIn(line)) {
                peptide = getSubUtilSimple(line, re_pep)
                protein = getSubUtilSimple(line, re_pro)
                transcript = protein.split("|").get(2)
                outLine = spectrum + "\t" + peptide + "\t" + protein + "\t" + transcript + "\n"
                appendFile(outLine, peptideTranscriptFilePath)
            }
            if (re_alter.toRegex().containsMatchIn(line)){
                protein = getSubUtilSimple(line,re_alter)
                transcript = protein.split("|").get(2)
                outLine = spectrum+"\t"+peptide+"\t"+protein+"\t"+transcript + "\n"
                appendFile(outLine,peptideTranscriptFilePath)
            }
        }
        return peptideTranscriptFilePath //返回文件名

    }
    fun generatePeptideCounts(): String {
        //对文件进行处理，得到相应肽段支持较多的转录本
        //val filename = ProjectPath.msDir+"peptideTranscript.txt"
        val filename = "/Users/toby/Desktop/neoantigenData/peptideTranscript_2.txt"
        var file = getFileLines(filename)
        val size = file.size
        //val peptideCountsFilePath = ProjectPath.msDir+"peptideCounts.txt"
        val peptideCountsFilePath = "/Users/toby/Desktop/neoantigenData/peptideCounts_2.txt"

        val spectrumCountsMap = HashMap<String,Double>()
        val transcriptCountsMap = HashMap<String,Double>()
        var outLine = ""

        /*
        //获取第四列transcript信息，并去重
        for (i in 0..size){
            var line = file.get(i)
            transcriptSets.add(line.split("\t").get(3))
        }

        //计算各个transcript counts信息
        while (transcriptSets_it.hasNext()) {
            var count = 0
            for (i in 0..size) {
                var line = file.get(i)
                if (transcriptSets_it.next().equals(line.split("\t").get(3))){
                    count ++
                }
            }
            transcriptCountsMap[transcriptSets_it.next()] = count
        }

        //return transcriptCountsMap
        //return？
        */

        for (i in 0..size-1){
            val line = file.get(i)
            when (spectrumCountsMap.containsKey(line.split("\t").get(0))){
                true -> spectrumCountsMap.replace(line.split("\t").get(0), spectrumCountsMap[line.split("\t").get(0)]!!.toDouble()+1)
                false -> spectrumCountsMap[line.split("\t").get(0)] = 1.0
            }
        }
        //直接使用map进行计算,没有考虑其他因素
        for (i in 0..size-1){
            var line = file.get(i)
            when (transcriptCountsMap.containsKey(line.split("\t").get(3))) {
                true -> transcriptCountsMap.replace(line.split("\t").get(3),transcriptCountsMap[line.split("\t").get(3)]!!.toDouble()+(1/spectrumCountsMap.getValue(line.split("\t").get(0))))
                false -> transcriptCountsMap[line.split("\t").get(3)] = 1/spectrumCountsMap.getValue(line.split("\t").get(0))
            }
        }

        transcriptCountsMap.entries.forEach {
            outLine = it.key+"\t"+it.value+"\n"
            appendFile(outLine,peptideCountsFilePath)
        }

        return peptideCountsFilePath

    }
//}