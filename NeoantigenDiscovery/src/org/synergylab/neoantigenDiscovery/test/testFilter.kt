package org.synergylab.neoantigenDiscovery.test

import org.synergylab.neoantigenDiscovery.utils.filterBindingAffinity
import org.synergylab.neoantigenDiscovery.utils.filterFPKM
import org.synergylab.neoantigenDiscovery.utils.filterGatkPass
import org.synergylab.neoantigenDiscovery.utils.filterPeptideCounts

fun main(args: Array<String>) {
    //filterGatkPass()
    //filterFPKM()
    //filterPeptideCounts()
    filterBindingAffinity()
}