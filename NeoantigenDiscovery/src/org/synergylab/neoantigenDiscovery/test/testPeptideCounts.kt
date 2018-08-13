package org.synergylab.neoantigenDiscovery.test

import org.synergylab.neoantigenDiscovery.NeoantigenMS.generatePeptideCounts
import org.synergylab.neoantigenDiscovery.NeoantigenMS.generatePeptideTranscriptFile

fun main(args: Array<String>) {
    //val test = generatePeptides
    //val test = peptideCounts
    generatePeptideTranscriptFile()
    val pepCounts = generatePeptideCounts()

}