package org.synergylab.neoantigenDiscovery.test

import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand
import org.synergylab.neoantigenDiscovery.utils.execute
import org.synergylab.neoantigenDiscovery.utils.text
import java.lang.Runtime
import java.util.Scanner

import java.io.IOException



   /* fun main(args: Array<String>) {
        try {
            val rec = RemoteExecuteCommand("192.168.6.203", "root", "rlibs402")
            rec.execute("perl /localdisk/software/scripts/exam.pl")
        } catch (e: IOException) {
            e.printStackTrace()
        }

        println("Finished")
    }
*/

    fun test(a:Int,b:Int): Int {
        val result = a+b
        return result

    }

    fun main(args: Array<String>) {
        //val rec = RemoteExecuteCommand("192.168.6.203", "root", "rlibs402")
        //执行命令
        //System.out.println(rec.execute("ifconfig"));
        //执行脚本
        //rec.execute("sh /usr/local/tomcat/bin/statup.sh");
        //rec.execute("perl /localdisk/software/scripts/exam.pl")
        //这个方法与上面最大的区别就是，上面的方法，不管执行成功与否都返回，
        //这个方法呢，如果命令或者脚本执行错误将返回空字符串
        //rec.executeSuccess("ifconfig");
        val command = "python2.7 /localdisk/software/seq2HLA2.2/seq2HLA.py -1 /localdisk/jobs/lung/180313002ML/cancer/cleandata/180313002ML_cancer_fq1_clean.fq " +
                "-2 /localdisk/jobs/lung/180313002ML/cancer/cleandata/180313002ML_cancer_fq2_clean.fq -r /localdisk/jobs/lung/180313002ML/cancer/hla/180313002ML -p 10"
        val proc = Runtime.getRuntime().exec(command)
        Scanner(proc.inputStream).use{
            while (it.hasNextLine()) println{it.nextLine()}
        }
        //val a = "hello"
        //val b = "${a.get(0)}_${a.get(1)}"
        //println(b)
        //val result = test(5,8)
        //println(result)
    }
