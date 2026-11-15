//
// MZMINE_PROCESS — stage 2: feature detection, alignment, gap filling, ion
// identity networking and exports, for all files of one polarity.
//
// Runs MZmine in batch mode with a shipped settings file
// (assets/mzmine/soil_lipidomics_{pos,neg}.mzbatch). Input files and output
// names are set on the command line (-i, -o), so the same settings file serves
// any dataset. MZmine 4 needs a user login: pass --mzmine_user with the path to
// an .mzuser file (created by logging in once, see docs/usage.md). The file is
// copied into a task-local home so nothing is written to the real home folder.
//
// Outputs are copied to fixed names so later stages do not depend on MZmine's
// export naming:
//   <prefix>.features_full.csv   all features, per-file area/height/mz/rt
//   <prefix>.gnps_quant.csv      GNPS FBMN quantification table   (if exported)
//   <prefix>.gnps.mgf            GNPS FBMN MS2 spectra             (if exported)
//   <prefix>.sirius.mgf          SIRIUS input                      (if exported)
//   <prefix>.annotations.csv     MZmine lipid annotations          (if exported)
//
process MZMINE_PROCESS {
    tag "${polarity}"
    label 'process_high'

    input:
    tuple val(polarity), path(mzmls, stageAs: 'mzml/*'), path(batch)

    output:
    tuple val(polarity), path("${prefix}.features_full.csv"), emit: features
    tuple val(polarity), path("${prefix}.gnps_quant.csv"),    emit: gnps_quant, optional: true
    tuple val(polarity), path("${prefix}.gnps.mgf"),          emit: gnps_mgf,   optional: true
    tuple val(polarity), path("${prefix}.sirius.mgf"),        emit: sirius_mgf, optional: true
    tuple val(polarity), path("${prefix}.annotations.csv"),   emit: annotations, optional: true
    path "${prefix}.mzmine.log",                               emit: log
    path 'versions.yml',                                       emit: versions

    script:
    prefix = "${params.run_name}_${polarity}"
    def mzmine = params.mzmine_cmd ?: 'mzmine'
    def user   = params.mzmine_user ? "-u \"\$PWD/mzmine_home/user.mzuser\"" : ''
    """
    ls -1 "\$PWD"/mzml/*.mzML > files.txt
    mkdir -p mzmine_out mzmine_tmp mzmine_home
    export JAVA_TOOL_OPTIONS="-Djdk.lang.Process.launchMechanism=FORK -Duser.home=\$PWD/mzmine_home \${JAVA_TOOL_OPTIONS:-}"
    ${params.mzmine_user ? "cp \"${params.mzmine_user}\" mzmine_home/user.mzuser" : ''}

    set +e
    ${mzmine} -b ${batch} -i files.txt -o "mzmine_out/${prefix}" -t mzmine_tmp --threads ${task.cpus} ${user} > ${prefix}.mzmine.log 2>&1
    status=\$?
    set -e
    tail -n 40 ${prefix}.mzmine.log
    if [ \$status -ne 0 ]; then
        if grep -q -E "needs to login again|Cannot read specified user file|no users logged" ${prefix}.mzmine.log; then
            echo "MZmine login missing or expired. Log in again once (mzmine --login-console)" >&2
            echo "and pass the new .mzuser file with --mzmine_user. See docs/usage.md (MZmine)." >&2
        fi
        echo "MZmine exited with status \$status; see ${prefix}.mzmine.log" >&2
        exit \$status
    fi

    copy_first() {
        found=\$(ls \$1 2>/dev/null | head -n 1)
        if [ -n "\$found" ]; then cp "\$found" "\$2"; fi
    }
    copy_first 'mzmine_out/*full*.csv'          ${prefix}.features_full.csv
    copy_first 'mzmine_out/*gnps*quant*.csv'    ${prefix}.gnps_quant.csv
    copy_first 'mzmine_out/*gnps*.mgf'          ${prefix}.gnps.mgf
    copy_first 'mzmine_out/*sirius*.mgf'        ${prefix}.sirius.mgf
    copy_first 'mzmine_out/*annotations*.csv'   ${prefix}.annotations.csv
    if [ ! -s ${prefix}.features_full.csv ]; then
        echo "MZmine finished but no full feature table was exported. Files written:" >&2
        ls -la mzmine_out >&2
        exit 1
    fi
    rm -rf mzmine_tmp mzmine_home

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mzmine: \$(grep -o -m1 'Starting mzmine [0-9.]*' ${prefix}.mzmine.log | sed 's/Starting mzmine //')
    END_VERSIONS
    """

    stub:
    prefix = "${params.run_name}_${polarity}"
    """
    echo "id,mz,rt" > ${prefix}.features_full.csv
    touch ${prefix}.mzmine.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mzmine: stub
    END_VERSIONS
    """
}
