package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess

import org.synergylab.neoantigenDiscovery.utils.*

//得到含有非同义突变的转录组
//得到方法还需确认

class MutatedTranscript(){
    //parameter

    //source data
    //val gmapFile = ProjectPath.cancerSampleDir+"gmap_out" //文件名暂定为gmap_out,地址还没定好
    //val NonSynonymousSomaticMutation = 0 //突变信息

    //class properties
    //val mutatedTranscript ＝ generateMutatedTranscript()


    //function
    fun generateGmapReference(gmapFile: String,outDir: String): String {
        var file = getFileLines(gmapFile)
        var size = file.size
        //val gmapReferenceFile = ProjectPath.cancerSampleDir+"gmapReference.txt"
        val gmapReferenceFile= "${outDir}gmapReference.txt"


        //gmap参数设置为-A
        val regexTrans = ">(.*?) len"
        val regexChr = "genome Chr(.*?):"
        val regexStrand = "\\((.*?) strand\\)"
        val regexPosition = "Chr\\d+:(\\d+)-(\\d+)  \\((\\d+)-(\\d+)\\)"
        val regexPath = "(Path .*?): query"

        var transcriptId: String = ""
        var path: String = ""
        var chr: String = ""
        var strand: String = ""
        var trinityStart: String = ""
        var trinityEnd: String = ""
        var refStart: String = ""
        var refEnd: String = ""
        var control: Boolean = false
        var outLine = ""

        for (i in 0 until size){
            val line = file.get(i)
            if (regexTrans.toRegex().containsMatchIn(line)){
                transcriptId = getSubUtilSimple(line, regexTrans)
            }
            if (regexPath.toRegex().containsMatchIn(line)){
                path = getSubUtilSimple(line,regexPath)
                control = when ("Path 1".equals(path)){
                    true -> true
                    false -> false
                }
            }
            if (control){
                if (regexChr.toRegex().containsMatchIn(line)){
                    chr = getSubUtilSimple(line,regexChr)
                }
                if (regexStrand.toRegex().containsMatchIn(line)){
                    strand = getSubUtilSimple(line,regexStrand)
                }
                if (regexPosition.toRegex().containsMatchIn(line)){
                    refStart = getSpecificSubUtilSimple(line,regexPosition,1)
                    refEnd = getSpecificSubUtilSimple(line,regexPosition,2)
                    trinityStart = getSpecificSubUtilSimple(line,regexPosition,3)
                    trinityEnd = getSpecificSubUtilSimple(line,regexPosition,4)
                    outLine = "${transcriptId}\t${chr}\t${strand}\t${trinityStart}\t${trinityEnd}\t${refStart}\t${refEnd}\n"
                    appendFile(outLine,gmapReferenceFile)
                }
            }
        }

        /*
        //gmap参数设为-Z
        val re_chr_info = "Chr(.*?):(\\d+)\\.\\.(\\d+)"
        val re_other_info = "(.*?):(\\d+)\\.\\.(\\d+)"
        var transcriptId: String = ""
        var info: String = ""
        var chr: String = ""
        var strand: String = ""
        var trans_start: String = ""
        var trans_end: String = ""
        var refStart: String = ""
        var refEnd: String = ""
        var outLine = ""

        for (i in 0..size-1){
            val line = file.get(i)
            val lineInfo = line.split(" ")
            if (line.startsWith(">")){
                transcriptId = lineInfo.get(0).substring(1)
                info = lineInfo.get(9)
                chr = getSpecificSubUtilSimple(info,re_chr_info,1)
                refStart = getSpecificSubUtilSimple(info,re_chr_info,2)
                refEnd = getSpecificSubUtilSimple(info,re_chr_info,3)
                strand = lineInfo.get(10)
                trans_start = lineInfo.get(7).split("..").get(0)
                trans_end = lineInfo.get(7).split("..").get(1)
                if (chr == ""){
                    chr = getSpecificSubUtilSimple(info,re_other_info,1)
                    refStart = getSpecificSubUtilSimple(info,re_other_info,2)
                    refEnd = getSpecificSubUtilSimple(info,re_other_info,3)
                }
                outLine = "${transcriptId}\t${chr}\t${strand}\t${trans_start}\t${tran_end}\t${refStart}\t${refEnd}\n"
                appendFile(outLine,gmapReferenceFile)
            }
        }
        */

        return gmapReferenceFile

    }

    fun generateMutatedTranscript(gmapReference: String,annotation: String,outDir: String): String {

        var gmapReferenceFile = getFileLines(gmapReference)
        val gmapReferenceSize = gmapReferenceFile.size
        //annovar文件要对突变类型进行筛选
        var annotationFile = getFileLines(annotation)
        val annotationSize = annotationFile.size
        //val MutatedTranscriptFile = ProjectPath.cancerSampleDir+"MutatedTranscript.txt"
        val MutatedTranscriptFile = "${outDir}MutatedTranscript.txt"
        val header = "TrinityID\tChr\tAlt_Base\tStrand\tTrinity_start\tTrinity_end\tRef_start\tRef_end\tPosition_start\tPosition_end\tTranscript\tMutationType\tProteinPosition\tAminoAcid\tCodons\n"
        appendFile(header, MutatedTranscriptFile)

        val regexLocation1 = "Chr(.*):(\\d+)-(\\d+)"
        val regexLocation2 = "Chr(.*):(\\d+)"
        var startAnno = 0
        var endAnno = 0
        var locationAnno = 0
        var outLine = ""


        loop@ for (i in 0 until gmapReferenceSize) {
            val gmapLine = gmapReferenceFile.get(i)
            val gmapLineInfo = gmapLine.split("\t")
            val transcriptID_gmap = gmapLineInfo.get(0)
            val chrGmap = gmapLineInfo.get(1)
            val strandGmap = gmapLineInfo.get(2)
            val trinityStartGmap = gmapLineInfo.get(3).toInt()
            val trinityEndGmap = gmapLineInfo.get(4).toInt()
            val refStartGmap = gmapLineInfo.get(5).toInt()
            val refEndGmap = gmapLineInfo.get(6).toInt()

            for (j in 0 until annotationSize) {
                val annotationLine = annotationFile.get(j)
                val annotationLineInfo = annotationLine.split("\t")
                if (!annotationLine.startsWith("#")) {
                    val chrom = annotationLineInfo.get(1)
                    val altBase = annotationLineInfo.get(2)
                    val chrAnno = getSpecificSubUtilSimple(chrom, regexLocation2, 1)
                    val transcript = annotationLineInfo.get(4)
                    val mutationType = annotationLineInfo.get(6)
                    val proteinPosition = annotationLineInfo.get(9)
                    val aminoAcid = annotationLineInfo.get(10)
                    val codons = annotationLineInfo.get(11)
                    val outInfo = "${transcript}\t${mutationType}\t${proteinPosition}\t${aminoAcid}\t${codons}"
                    if (chrGmap.equals(chrAnno)) {
                        if (regexLocation1.toRegex().containsMatchIn(chrom)) {
                            startAnno = getSpecificSubUtilSimple(chrom, regexLocation1, 2).toInt()
                            endAnno = getSpecificSubUtilSimple(chrom, regexLocation1, 3).toInt()
                            if ("+".equals(strandGmap) && (refStartGmap <= startAnno) && (endAnno <= refEndGmap)) {
                                outLine = "${transcriptID_gmap}\t${chrGmap}\t${altBase}\t${strandGmap}\t${trinityStartGmap}\t${trinityEndGmap}\t${refStartGmap}\t${refEndGmap}\t${startAnno}\t${endAnno}\t${outInfo}\n"
                                appendFile(outLine, MutatedTranscriptFile)
                                continue@loop
                            } else if ("-".equals(strandGmap) && (refEndGmap <= startAnno) && (endAnno <= refStartGmap)) {
                                outLine = "${transcriptID_gmap}\t${chrGmap}\t${altBase}\t${strandGmap}\t${trinityStartGmap}\t${trinityEndGmap}\t${refStartGmap}\t${refEndGmap}\t${startAnno}\t${endAnno}\t${outInfo}\n"
                                appendFile(outLine, MutatedTranscriptFile)
                                continue@loop
                            }
                        } else if (regexLocation2.toRegex().containsMatchIn(chrom)) {
                            locationAnno = getSpecificSubUtilSimple(chrom, regexLocation2, 2).toInt()
                            startAnno = locationAnno
                            val endReplaceAnno = "-"
                            if ("+".equals(strandGmap) && (refStartGmap <= locationAnno) && (locationAnno <= refEndGmap)) {
                                outLine = "${transcriptID_gmap}\t${chrGmap}\t${altBase}\t${strandGmap}\t${trinityStartGmap}\t${trinityEndGmap}\t${refStartGmap}\t${refEndGmap}\t${startAnno}\t${endAnno}\t${outInfo}\n"
                                appendFile(outLine, MutatedTranscriptFile)
                                continue@loop
                            } else if ("-".equals(strandGmap) && (refEndGmap <= locationAnno) && (locationAnno <= refStartGmap)) {
                                outLine = "${transcriptID_gmap}\t${chrGmap}\t${altBase}\t${strandGmap}\t${trinityStartGmap}\t${trinityEndGmap}\t${refStartGmap}\t${refEndGmap}\t${startAnno}\t${endAnno}\t${outInfo}\n"
                                appendFile(outLine, MutatedTranscriptFile)
                                continue@loop
                            }
                        }
                    }
                }
            }
        }
        return MutatedTranscriptFile
    }

}