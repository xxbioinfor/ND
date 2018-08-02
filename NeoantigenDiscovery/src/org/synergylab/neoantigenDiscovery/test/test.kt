package org.synergylab.neoantigenDiscovery.test

import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand
import org.synergylab.neoantigenDiscovery.utils.execute
import org.synergylab.neoantigenDiscovery.utils.text

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

    fun main(args: Array<String>) {
        val rec = RemoteExecuteCommand("192.168.6.203", "root", "rlibs402")
        //执行命令
        //System.out.println(rec.execute("ifconfig"));
        //执行脚本
        //rec.execute("sh /usr/local/tomcat/bin/statup.sh");
        rec.execute("perl /localdisk/software/scripts/exam.pl")
        //这个方法与上面最大的区别就是，上面的方法，不管执行成功与否都返回，
        //这个方法呢，如果命令或者脚本执行错误将返回空字符串
        //rec.executeSuccess("ifconfig");

    }
