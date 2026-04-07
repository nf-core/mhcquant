process PRIDE_FETCH_SDRF {
    label 'process_single'
    tag "$pride_id"

    conda "bioconda::pridepy=0.0.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pridepy:0.0.12--pyhdfd78af_0' :
        'quay.io/biocontainers/pridepy:0.0.12--pyhdfd78af_0' }"

    input:
    val pride_id

    output:
    path "*.sdrf.tsv" , emit: sdrf
    tuple val("${task.process}"), val('pridepy'), eval("pip show pridepy 2>/dev/null | grep Version | cut -d' ' -f2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat <<'EOF' > fetch_sdrf.py
import json
import subprocess
import sys

pride_id = "${pride_id}"

subprocess.run(
    ["pridepy", "stream-files-metadata", "-a", pride_id, "-o", "files_metadata.json"],
    check=True
)

with open("files_metadata.json") as f:
    files = json.load(f)

sdrf_files = [f for f in files if f["fileName"].endswith(".sdrf.tsv")]
if not sdrf_files:
    print(f"ERROR: No SDRF file found for PRIDE project {pride_id}", file=sys.stderr)
    print("Please provide a local SDRF file using --input <path_to_sdrf.tsv>", file=sys.stderr)
    sys.exit(1)

if len(sdrf_files) > 1:
    names = [f["fileName"] for f in sdrf_files]
    print(f"WARNING: Multiple SDRF files found: {names}. Using {names[0]}", file=sys.stderr)

sdrf_name = sdrf_files[0]["fileName"]
print(f"Found SDRF file: {sdrf_name}")
subprocess.run(
    ["pridepy", "download-file-by-name", "-a", pride_id, "-f", sdrf_name, "-o", ".", "-p", "ftp"],
    check=True
)
EOF

    python3 fetch_sdrf.py
    """

    stub:
    """
    touch ${pride_id}.sdrf.tsv
    """
}
