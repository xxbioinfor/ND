package org.synergylab.neoantigenDiscovery.utils

import java.io.BufferedReader
import java.io.InputStreamReader

fun String.execute(): Process{
    val runtime = Runtime.getRuntime()
    return runtime.exec(this)
}

fun Process.text(): String{
    var output=""
    val inputStream = this.inputStream
    val isr = InputStreamReader(inputStream)
    val reader = BufferedReader(isr)
    var line: String? = ""
    while(line != null){
        line = reader.readLine()
        output += line + "\n"
    }
    return output
}

