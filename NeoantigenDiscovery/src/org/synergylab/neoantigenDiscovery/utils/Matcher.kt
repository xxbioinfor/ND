package org.synergylab.neoantigenDiscovery.utils

import java.util.regex.Pattern
import com.sun.tools.internal.xjc.model.Multiplicity.group
import java.util.ArrayList


/**
 * 正则表达式匹配两个指定字符串中间的内容
 * @param source
 * @return
 */
fun getSubUtil(source: String, rgex: String): List<String> {
    val list = ArrayList<String>()
    val pattern = Pattern.compile(rgex)// 匹配的模式
    val m = pattern.matcher(source)
    while (m.find()) {
        var i = 1
        list.add(m.group(i))
        i++
    }
    return list
}

/**
 * 返回单个字符串，若匹配到多个的话就返回第一个，方法与getSubUtil一样
 * @param source
 * @param rgex
 * @return
 */
fun getSubUtilSimple(source: String, rgex: String): String{
    val pattern = Pattern.compile(rgex)// 匹配的模式
    val m = pattern.matcher(source)
    while (m.find()) {
        return m.group(1)
    }
    return ""
}

/**
 * 返回单个字符串，并根据所匹配的位置，输出相应的信息
 * @param source
 * @param rgex
 * @param num
 * @return
 */
fun getSpecificSubUtilSimple(source: String, rgex: String, num: Int): String {
    val pattern = Pattern.compile(rgex)// 匹配的模式
    val m = pattern.matcher(source)
    while (m.find()) {
        return m.group(num)
    }
    return ""
}



class Matcher(source: String){
    val sourceString = source

    var sampleID = sampleIDMatcher()

    /**
     * 由原始测序数据生成样本ID，格式为测序时间序列号＋癌症类型
     * @author wxx
     * @since  V0.1
     */
    private fun sampleIDMatcher(): String {

        val pattern = """(\d*?)\D.*"""
        val matcher = Pattern.compile(pattern).matcher(this.sourceString)

        while (matcher.find()){
            sampleID =  matcher.group(1)+ProjectPath.cancerType
    }
        return sampleID
    }


}