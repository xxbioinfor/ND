package NeoantigenDiscovery

//存储一些data class
data class NeoantigenResult(
        val id:Int,
        val peptide: String,
        val HLAtype: String,
        val Affinity: Double,
        val rank: Double
)

data class BindingResult(
        val posID: Int,
        val peptide: String,
        val HLAtype: String,
        val Affinity: Double,
        val rank: Double,
        val core: String,
        val NumBinders: Int
)

data class PeptideTranscriptResult(
        val spectrum: String,
        val peptide: String,
        val protein: String,
        val transcript: String
)