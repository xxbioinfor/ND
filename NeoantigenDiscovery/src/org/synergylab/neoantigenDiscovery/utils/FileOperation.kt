package org.synergylab.neoantigenDiscovery.utils

import java.io.File
import java.nio.charset.Charset

fun getFileLines(filename: String): List<String> {
    return File(filename).readLines(Charset.forName("UTF-8"))
}

fun getFileContent(filename: String): String {
    val f = File(filename)
    return f.readText(Charset.forName("UTF-8"))
}

fun writeFile(text: String, destFile: String) {
    val f = File(destFile)
    if (!f.exists()) {
        f.createNewFile()
    }
    f.writeText(text, Charset.defaultCharset())
}

fun appendFile(text: String, destFile: String) {
    val f = File(destFile)
    if (!f.exists()) {
        f.createNewFile()
    }
    f.appendText(text, Charset.defaultCharset())
}
