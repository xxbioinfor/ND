package org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess

import org.synergylab.neoantigenDiscovery.utils.*

//得到含有非同义突变的转录组
//得到方法还需确认

//class MutatedTranscript(){
    //parameter

    //source data
    //val GmapReference = 0 //文件名暂定为gmap_out
    //val NonSynonymousSomaticMutation = 0 //突变信息

    //class properties
    //val mutatedTranscript ＝ generateMutatedTranscript()


    //function
    fun generateGmapReference(): String {
        //val filename = ProjectPath.cancerSampleDir+"gmap_out" //地址还没定好
        val filename = "/Users/toby/Desktop/neoantigenData/gmap_out"
        var file = getFileLines(filename)
        var size = file.size
        //val gmapReferenceFile = ProjectPath.cancerSampleDir+"gmapReference.txt"
        val gmapReferenceFile= "/Users/toby/Desktop/neoantigenData/gmapReference_all2.txt"

        /*
        //gmap参数设置为-A
        val re_trans = ">(.*?) len"
        val re_chr = "genome Chr(.*?):"
        val re_strand = "\\((.*?) strand\\)"
        val re_ref = "Chr\\d+:(\\d+)-(\\d+)"
        val re_triniy = "\\((\\d+)-(\\d+)\\)"
        val re_path = "(Path .*?): query"

        var transcriptId: String = ""
        var path: String = ""
        var chr: String = ""
        var strand: String = ""
        var trans_start: String = ""
        var trans_end: String = ""
        var ref_start: String = ""
        var ref_end: String = ""
        var control: Boolean = false
        var outLine = ""



        for (i in 0..size-1){
            val line = file.get(i)
            if (re_trans.toRegex().containsMatchIn(line)){
                transcriptId = getSubUtilSimple(line, re_trans)
            }
            if (re_path.toRegex().containsMatchIn(line)){
                path = getSubUtilSimple(line,re_path)
                control = when (path.equals("Path 1")){
                    true -> true
                    false -> false
                }
            }
            if (control){
                if (re_chr.toRegex().containsMatchIn(line)){
                    chr = getSubUtilSimple(line,re_chr)
                }
                if (re_strand.toRegex().containsMatchIn(line)){
                    strand = getSubUtilSimple(line,re_strand)
                }
                if (re_ref.toRegex().containsMatchIn(line)){
                    ref_start = getSpecificSubUtilSimple(line,re_ref,1)
                    ref_end = getSpecificSubUtilSimple(line,re_ref,2)
                }
                if (re_triniy.toRegex().containsMatchIn(line)){
                    trans_start = getSpecificSubUtilSimple(line,re_triniy,1)
                    trans_end = getSpecificSubUtilSimple(line,re_triniy,2)
                    outLine = transcriptId+"\t"+chr+"\t"+strand+"\t"+trans_start+"\t"+trans_end+"\t"+ref_start+"\t"+ref_end+"\n"
                    appendFile(outLine,gmapReferenceFile)
                    control = false
                }
            }
        }
        */

        //gmap参数设为-Z
        val re_chr_info = "Chr(.*?):(\\d+)\\.\\.(\\d+)"
        val re_other_info = "(.*?):(\\d+)\\.\\.(\\d+)"
        var transcriptId: String = ""
        var info: String = ""
        var chr: String = ""
        var strand: String = ""
        var trans_start: String = ""
        var trans_end: String = ""
        var ref_start: String = ""
        var ref_end: String = ""
        var outLine = ""

        for (i in 0..size-1){
            val line = file.get(i)
            if (line.startsWith(">")){
                transcriptId = line.split(" ").get(0).substring(1)
                info = line.split(" ").get(9)
                chr = getSpecificSubUtilSimple(info,re_chr_info,1)
                ref_start = getSpecificSubUtilSimple(info,re_chr_info,2)
                ref_end = getSpecificSubUtilSimple(info,re_chr_info,3)
                strand = line.split(" ").get(10)
                trans_start = line.split(" ").get(7).split("..").get(0)
                trans_end = line.split(" ").get(7).split("..").get(1)
                if (chr == ""){
                    chr = getSpecificSubUtilSimple(info,re_other_info,1)
                    ref_start = getSpecificSubUtilSimple(info,re_other_info,2)
                    ref_end = getSpecificSubUtilSimple(info,re_other_info,3)
                }
                outLine = transcriptId+"\t"+chr+"\t"+strand+"\t"+trans_start+"\t"+trans_end+"\t"+ref_start+"\t"+ref_end+"\n"
                appendFile(outLine,gmapReferenceFile)
            }
        }

        return gmapReferenceFile

    }

    fun generateMutatedTranscript(): String {
        //val gmapFileName = ProjectPath.cancerSampleDir+"gmapReference.txt"
        val gmapFileName = "/Users/toby/Desktop/neoantigenData/gmapReference_all2.txt"
        var gmapFile = getFileLines(gmapFileName)
        val gmapSize = gmapFile.size
        //annovar文件要对突变类型进行筛选
        //val annotationFileName = ProjectPath.sampleBaseDir+"somatic_mutations_pass.txt"
        val annotationFileName = "/Users/toby/Desktop/neoantigenData/filter_pass_vep.vcf"
        var annotationFile = getFileLines(annotationFileName)
        val annotationSize = annotationFile.size
        //val MutatedTranscriptFile = ProjectPath.cancerSampleDir+"MutatedTranscript.txt"
        val MutatedTranscriptFile = "/Users/toby/Desktop/neoantigenData/MutatedTranscript_all.txt"
        val header = "TrinityID\tChr\tAlt_Base\tStrand\tTrinity_start\tTrinity_end\tRef_start\tRef_end\tPosition_start\tPosition_end\tTranscript\tMutationType\tProteinPosition\tAminoAcid\tCodons\n"
        appendFile(header, MutatedTranscriptFile)

        val re_location_1 = "Chr(.*):(\\d+)-(\\d+)"
        val re_location_2 = "Chr(.*):(\\d+)"
        var chr_anno = ""
        var start_anno = 0
        var end_anno = 0
        var location_anno = 0
        var outLine = ""


        loop@ for (i in 0..gmapSize - 1) {
            var gmapLine = gmapFile.get(i)
            var transcriptID_gmap = gmapLine.split("\t").get(0)
            var chr_gmap = gmapLine.split("\t").get(1)
            var strand_gmap = gmapLine.split("\t").get(2)
            var trinity_start_gmap = gmapLine.split("\t").get(3).toInt()
            var trinity_end_gmap = gmapLine.split("\t").get(4).toInt()
            var ref_start_gmap = gmapLine.split("\t").get(5).toInt()
            var ref_end_gmap = gmapLine.split("\t").get(6).toInt()

            for (j in 41..annotationSize - 1) {
                var annotationLine = annotationFile.get(j)
                if (!annotationLine.startsWith("#")) {
                    val chrom = annotationLine.split("\t").get(1)
                    val alt_base = annotationLine.split("\t").get(2)
                    val chr_anno = getSpecificSubUtilSimple(chrom, re_location_2, 1)
                    val transcript = annotationLine.split("\t").get(4)
                    val mutationType = annotationLine.split("\t").get(6)
                    val proteinPosition = annotationLine.split("\t").get(9)
                    val aminoAcid = annotationLine.split("\t").get(10)
                    val codons = annotationLine.split("\t").get(11)
                    val outInfo = transcript+"\t"+mutationType+"\t"+proteinPosition+"\t"+aminoAcid+"\t"+codons
                    if (chr_gmap.equals(chr_anno)) {
                        if (re_location_1.toRegex().containsMatchIn(chrom)) {
                            start_anno = getSpecificSubUtilSimple(chrom, re_location_1, 2).toInt()
                            end_anno = getSpecificSubUtilSimple(chrom, re_location_1, 3).toInt()
                            if (strand_gmap.equals("+") && (ref_start_gmap <= start_anno) && (end_anno <= ref_end_gmap)) {
                                outLine = transcriptID_gmap+"\t"+chr_gmap+"\t"+alt_base+"\t"+strand_gmap+"\t"+trinity_start_gmap+"\t"+trinity_end_gmap+"\t"+ref_start_gmap+"\t"+ref_end_gmap+"\t"+start_anno+"\t"+end_anno+"\t"+outInfo+"\n"
                                appendFile(outLine, MutatedTranscriptFile)
                                continue@loop
                            } else if (strand_gmap.equals("-") && (ref_end_gmap <= start_anno) && (end_anno <= ref_start_gmap)) {
                                outLine = transcriptID_gmap+"\t"+chr_gmap+"\t"+alt_base+"\t"+strand_gmap+"\t"+trinity_start_gmap+"\t"+trinity_end_gmap+"\t"+ref_start_gmap+"\t"+ref_end_gmap+"\t"+start_anno+"\t"+end_anno+"\t"+outInfo+"\n"
                                appendFile(outLine, MutatedTranscriptFile)
                                continue@loop
                            }
                        } else if (re_location_2.toRegex().containsMatchIn(chrom)) {
                            location_anno = getSpecificSubUtilSimple(chrom, re_location_2, 2).toInt()
                            start_anno = location_anno
                            val end_replace_anno = "-"
                            if (strand_gmap.equals("+") && (ref_start_gmap <= location_anno) && (location_anno <= ref_end_gmap)) {
                                outLine = transcriptID_gmap+"\t"+chr_gmap+"\t"+alt_base+"\t"+strand_gmap+"\t"+trinity_start_gmap+"\t"+trinity_end_gmap+"\t"+ref_start_gmap+"\t"+ref_end_gmap+"\t"+start_anno+"\t"+end_replace_anno+"\t"+outInfo+"\n"
                                appendFile(outLine, MutatedTranscriptFile)
                                continue@loop
                            } else if (strand_gmap.equals("-") && (ref_end_gmap <= location_anno) && (location_anno <= ref_start_gmap)) {
                                outLine = transcriptID_gmap+"\t"+chr_gmap+"\t"+alt_base+"\t"+strand_gmap+"\t"+trinity_start_gmap+"\t"+trinity_end_gmap+"\t"+ref_start_gmap+"\t"+ref_end_gmap+"\t"+start_anno+"\t"+end_replace_anno+"\t"+outInfo+"\n"
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

//}