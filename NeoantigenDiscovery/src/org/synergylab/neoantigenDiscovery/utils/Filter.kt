package org.synergylab.neoantigenDiscovery.utils

//进行一些过滤操作
//涉及一些过滤参数
//

//fun transcript expression
//fun

class Filter() {
    //parameters

    //source data
    //val mutect2File = ProjectPath.sampleBaseDir+"somatic_mutations.vcf"
    //val mutect2File = "/Users/toby/Desktop/neoantigenData/somatic_mutations.vcf"
    //class properties

    //function
    fun filterGatkPass(mutationFile: String,outDir: String):String {

        val file = getFileLines(mutationFile)
        val size = file.size
        val passFile = outDir+"somatic_mutations_pass.txt"
        val vepSourceFile = outDir+"somatic_mutations_pass.vcf"

        for (i in 0..size-1) {
            val line = file.get(i)
            var outLine = ""
            if (line.startsWith("#")) {
                outLine = line + "\n"
                appendFile(outLine, vepSourceFile)
            }
            else if (line.split("\t").get(6).equals("PASS")) {
                outLine = line + "\n"
                //val outLine = line.split("\t").get(0)+"\t"+line.split("\t").get(1)+"\t"+line.split("\t").get(2)+"\t"+line.split("\t").get(3)+"\t"+line.split("\t").get(4)+"\t"+line.split("\t").get(5)+"\t"+line.split("\t").get(6)+"\t"+line.split("\t").get(7)+"\t"+line.split("\t").get(8)+"\t"+line.split("\t").get(9)+"\n"
                appendFile(outLine, passFile)
                appendFile(outLine, vepSourceFile)
            }
        }
        return vepSourceFile
    }

    fun filterRsemFpkm(transcriptFile: String, rsemFile: String, minFpkm: Double, outDir: String): String {
        //val rsemFileName = ProjectPath.cancerSampleDir+"trinity_out_dir/RSEM.isoforms.results"
        val rsemFileName = rsemFile
        val trinityFileName = transcriptFile
        val minFpkm = minFpkm

        val rsemFile = getFileLines(rsemFileName)
        val rsemSize = rsemFile.size
        //val trinityFileName = ProjectPath.cancerSampleDir+"trinity_out_dir/Trinity.fasta"
        val trinityFile = getFileLines(trinityFileName)
        val trinitySize = trinityFile.size
        //val filterFile = ProjectPath.cancerSampleDir+"trinity_out_dir/Trinity_FPKM1.0.fasta"
        val filterFile = transcriptFile.split(".").get(0)+"_FPKM.txt"

        loop@ for (i in 0 until trinitySize) {
            val trinityLine = trinityFile.get(i)
            val trinityID = trinityLine.split("\t").get(0)

            for (j in 0 until rsemSize) {
                val rsemLine = rsemFile.get(j)
                val rsemID = rsemLine.split("\t").get(0)
                val fpkm = rsemLine.split("\t").get(6).toDouble()
                if (trinityID.equals(rsemID) && fpkm >= minFpkm) {
                    appendFile(trinityLine + "\t" + fpkm + "\n", filterFile)
                    continue@loop
                }
            }
        }
        return filterFile

        /*
        //使用trinity.fasta直接进行filter
        val fastaMap = HashMap<String, String>()
        var header = ""
        var fasta = ""

        for (i in 0 until trinitySize) {
            val trinityLine = trinityFile.get(i)
            if (trinityLine.startsWith(">")) {
                header = getSubUtilSimple(trinityLine, ">(.*?) len")
                fasta = ""
            } else {
                fasta += trinityLine
                fastaMap[header] = fasta
            }
        }

        loop@ for (i in 0 until trinitySize) {
            val trinityLine = trinityFile.get(i)
            if (trinityLine.startsWith(">")) {
                val trinityID = getSubUtilSimple(trinityLine, ">(.*?) len")
                for (j in 0 until rsemSize) {
                    val rsemLine = rsemFile.get(j)
                    val rsemID = rsemLine.split("\t").get(0)
                    val fpkm = rsemLine.split("\t").get(6).toDouble()
                    if (trinityID.equals(rsemID) && fpkm >= minFpkm) {
                        appendFile(trinityLine + "\n", filterFile)
                        appendFile(fastaMap[rsemID].toString() + "\n", filterFile)
                        continue@loop
                    }
                }
            }
        }
        */

    }

    fun filterPeptideCounts(transcriptFile: String, peptideCountsFile: String, minCounts: Double, outDir: String): String {
        //val peptideCountsFileName = ProjectPath.msDir+"peptideCounts.txt"
        val peptideCountsFileName = peptideCountsFile
        val mutatedTranscriptFileName = transcriptFile
        val minCounts = minCounts
        //val peptideCountsFileName = "/Users/toby/Desktop/neoantigenData/peptideCounts_2.txt"
        val peptideCountsFile = getFileLines(peptideCountsFileName)
        val peptideCountsSize = peptideCountsFile.size
        //val mutatedTranscriptFileName = ProjectPath.cancerSampleDir+"MutatedTranscript.txt"
        //val mutatedTranscriptFileName = "/Users/toby/Desktop/neoantigenData/MutatedTranscript_all.txt"
        val mutatedTranscriptFile = getFileLines(mutatedTranscriptFileName)
        val mutatedTranscriptSize = mutatedTranscriptFile.size

        val filterFile = transcriptFile.split(".").get(0)+"_pepCounts.txt"

        loop@ for (i in 0 until mutatedTranscriptSize) {
            val mutatedTranscriptLine = mutatedTranscriptFile.get(i)
            val trinityID = mutatedTranscriptLine.split("\t").get(0)

            for (j in 0 until peptideCountsSize) {
                val peptideLine = peptideCountsFile.get(j)
                val peptideID = peptideLine.split("\t").get(0)
                val peptideCounts = peptideLine.split("\t").get(1).toDouble()
                if (trinityID.equals(peptideID) && peptideCounts > minCounts) {
                    appendFile(mutatedTranscriptLine + "\t" + peptideCounts + "\n", filterFile)
                    continue@loop
                }
            }
        }

        return filterFile
    }

    fun filterBindingAffinity(affinityFile: String,hlaType: String,outDir: String): String {
        //file,hlaType,peptideLength
        val bindingAffinityFile = getFileLines(affinityFile)
        val bindingAffinitySize = bindingAffinityFile.size

        var count = hlaType.count { it == ',' } +1
        val bindingAffinityMap = HashMap<String, String>()

        val filterOutFile = outDir+"180313002ML_binding_filter.txt"

        if (count == 1) {
            for (i in 2 until bindingAffinitySize) {
                val bindingAffinityLine = bindingAffinityFile.get(i)
                val pos = bindingAffinityLine.split("\t").get(0)
                val peptide = bindingAffinityLine.split("\t").get(1)
                val id = bindingAffinityLine.split("\t").get(2)
                val newID = id.split("_").get(0) + "_" + id.split("_").get(1) + "_" + pos
                val bindingScore = bindingAffinityLine.split("\t").get(3)
                bindingAffinityMap[newID] = peptide + "\t" + bindingScore
            }

            for (i in 2 until bindingAffinitySize) {
                val bindingAffinityLine = bindingAffinityFile.get(i)
                val pos = bindingAffinityLine.split("\t").get(0)
                val peptide = bindingAffinityLine.split("\t").get(1)
                val id = bindingAffinityLine.split("\t").get(2)
                val idNumber = id.split("_").get(1)
                val newID = id.split("_").get(0) + "_" + id.split("_").get(1) + "_" + pos
                val bindingScore = bindingAffinityLine.split("\t").get(3).toDouble()

                if (newID.startsWith("mutated") && bindingScore < 500.0) {
                    val wildID = "wild_" + idNumber + "_" + pos
                    val wildBindingAffinity = bindingAffinityMap[wildID]!!.split("\t").get(1).toDouble()
                    val foldChange = wildBindingAffinity/bindingScore
                    if (wildBindingAffinity >= 500.0) {
                        val wildPeptide = bindingAffinityMap[wildID]!!.split("\t").get(0)
                        val outLine = newID + "\t" + peptide + "\t" + bindingScore + "\t" + wildID + "\t" + wildPeptide + "\t" + wildBindingAffinity + "\t" + foldChange + "\n"
                        appendFile(outLine, filterOutFile)
                    }
                }
            }
        }
        return filterOutFile
    }
}