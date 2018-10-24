package org.synergylab.neoantigenDiscovery.test

import org.synergylab.neoantigenDiscovery.NeoantigenTranscript.PostProcess.MutatedTranscript


fun main(args: Array<String>) {
    //println(" Generate Gmap Reference... ")
    val gmapReference = MutatedTranscript().generateGmapReference("/Users/toby/Desktop/neoantigenData/gmap_test_a")
    //val mutatedTranscript = MutatedTranscript().generateMutatedTranscript(gmapReference)

}