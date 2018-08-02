package org.synergylab.neoantigenDiscovery.NeoantigenMS

import org.synergylab.neoantigenDiscovery.utils.ProjectPath
import org.synergylab.neoantigenDiscovery.utils.RemoteExecuteCommand

//将质谱raw数据进行格式转换
//使用msConvert将raw转换成mzXML文件
//使用x!tandem将mzXML转换成tandem文件


//class MSFile(){
    //parameter

    //source data
    //val rawFile = 0 //来自cleanfileprocess处理得到的bamfile

    //class properties
    //val convertFile ＝ generateConvertFile()

    //function
    fun generateMzxmlFile(){
        val rec = RemoteExecuteCommand(ProjectPath.winIP,ProjectPath.winName, ProjectPath.winPwd)
        rec.execute("msconvert *.raw --mzXML -o dir")//需要修改
        rec.execute("scp *.mzXML root@"+ProjectPath.linuxIP+":"+ProjectPath.msDir)
        rec.execute("rlibs402")

        val re = RemoteExecuteCommand(ProjectPath.linuxIP,ProjectPath.linuxName,ProjectPath.linuxPwd)
        re.execute("xtandem")
        //生成tandem.pep.xml文件

    }

//}
