package org.synergylab.neoantigenDiscovery.utils

//进行一些过滤操作
//涉及一些过滤参数
//

//fun transcript expression
//fun

fun filterGatkPass(){
    //val fileName = ProjectPath.sampleBaseDir+"somatic_mutations.vcf"
    val fileName = "/Users/toby/Desktop/neoantigenData/somatic_mutations.vcf"
    val file = getFileLines(fileName)
    val size = file.size
    //val passFile = ProjectPath.sampleBaseDir+"somatic_mutations_pass.txt"
    //val vepSourceFile = ProjectPath.sampleBaseDir+"somatic_mutations_pass_vepsource.vcf"
    val passFile = "/Users/toby/Desktop/neoantigenData/somatic_mutations_pass.txt"
    val vepSourceFile = "/Users/toby/Desktop/neoantigenData/somatic_mutations_pass_vepsource.vcf"

    for (i in 0..size-1){
        val line = file.get(i)
        if (line.startsWith("#")){
            val outLine = line + "\n"
            appendFile(outLine,vepSourceFile)
        }
        if (!line.startsWith("#") && line.split("\t").get(6).equals("PASS")){
            val outline = line + "\n"
            //val outLine = line.split("\t").get(0)+"\t"+line.split("\t").get(1)+"\t"+line.split("\t").get(2)+"\t"+line.split("\t").get(3)+"\t"+line.split("\t").get(4)+"\t"+line.split("\t").get(5)+"\t"+line.split("\t").get(6)+"\t"+line.split("\t").get(7)+"\t"+line.split("\t").get(8)+"\t"+line.split("\t").get(9)+"\n"
            appendFile(outline,passFile)
            appendFile(outline,vepSourceFile)
        }
    }
}

fun filterFPKM(){
    //val rsemFileName = ProjectPath.cancerSampleDir+"trinity_out_dir/RSEM.isoforms.results"
    val rsemFileName = "/Users/toby/Desktop/neoantigenData/RSEM.isoforms.results"
    val rsemFile = getFileLines(rsemFileName)
    val rsemSize = rsemFile.size
    //val trinityFileName = ProjectPath.cancerSampleDir+"trinity_out_dir/Trinity.fasta"
    val trinityFileName = "/Users/toby/Desktop/neoantigenData/Trinity.fasta"
    val trinityFile = getFileLines(trinityFileName)
    val trinitySize = trinityFile.size
    //val filterFile = ProjectPath.cancerSampleDir+"trinity_out_dir/Trinity_FPKM1.0.fasta"
    val filterFile = "/Users/toby/Desktop/neoantigenData/Trinity_FPKM1.0.fasta"

    val fastaMap = HashMap<String,String>()
    var header = ""
    var fasta = ""

    for (i in 0..trinitySize-1){
        val trinityLine = trinityFile.get(i)
        if (trinityLine.startsWith(">")){
            header = getSubUtilSimple(trinityLine,">(.*?) len")
            fasta = ""
        }else{
            fasta += trinityLine
            fastaMap[header] = fasta
        }
    }

    loop@ for (i in 0..trinitySize-1){
        val trinityLine = trinityFile.get(i)
        if (trinityLine.startsWith(">")) {
            val trinityID = getSubUtilSimple(trinityLine, ">(.*?) len")
            for (j in 0..rsemSize - 1) {
                val rsemLine = rsemFile.get(j)
                val rsemID = rsemLine.split("\t").get(0)
                val fpkm = rsemLine.split("\t").get(6)
                if (trinityID.equals(rsemID) && fpkm > "1") {
                    appendFile(trinityLine+"\n", filterFile)
                    appendFile(fastaMap[rsemID].toString() + "\n", filterFile)
                    continue@loop
                }
            }
        }
    }

}

fun filterPeptideCounts(){
    //val peptideCountsFileName = ProjectPath.msDir+"peptideCounts.txt"
    val peptideCountsFileName = "/Users/toby/Desktop/neoantigenData/peptideCounts_2.txt"
    val peptideCountsFile = getFileLines(peptideCountsFileName)
    val peptideCountsSize = peptideCountsFile.size
    //val mutatedTranscriptFileName = ProjectPath.cancerSampleDir+"MutatedTranscript.txt"
    val mutatedTranscriptFileName = "/Users/toby/Desktop/neoantigenData/MutatedTranscript_all.txt"
    val mutatedTranscriptFile = getFileLines(mutatedTranscriptFileName)
    val mutatedTranscriptSize = mutatedTranscriptFile.size

    val filterFile = "/Users/toby/Desktop/neoantigenData/MutatedTranscript_pepCounts1.0.txt"

    loop@ for (i in 0..mutatedTranscriptSize-1){
        val mutatedTranscriptLine = mutatedTranscriptFile.get(i)
        val trinityID = mutatedTranscriptLine.split("\t").get(0)

        for (j in 0..peptideCountsSize-1){
            val peptideLine = peptideCountsFile.get(j)
            val peptideID = peptideLine.split("\t").get(0)
            val peptideCounts = peptideLine.split("\t").get(1)
            if (trinityID.equals(peptideID) && peptideCounts>"1"){
                appendFile(mutatedTranscriptLine+"\n",filterFile)
                continue@loop
            }
        }
    }
}