//
// libSignature() — a short hash of the vendored Python package lib/soillifeatlas.
//
// Nextflow's task hash covers the task script and the bin/ scripts it calls,
// but not code those scripts import. Modules that run the package put this
// hash in their script (and in versions.yml), so editing the package
// invalidates their cached results under -resume and every result records
// which package code produced it.
//
def libSignature() {
    def root = "${projectDir}/lib/soillifeatlas"
    def text = files("${root}/**.py")
        .sort { it.toString() }
        .collect { f -> f.toString().substring(root.length()) + '\n' + f.text }
        .join('\n')
    return text.md5().take(12)
}
