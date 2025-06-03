process EPICORE {
	tag "$meta.id"
	label 'process_high'

	conda "bioconda::epicore=0.1.3"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
		'https://depot.galaxyproject.org/singularity/epicore:10.1.3--pyhdfd78af_0' :
		'biocontainers/epicore:0.1.3--pyhdfd78af_0' }"

	input: 
		tuple val(meta_fasta), path(fasta)
		path(quantification_tsv)
		path(general_stats) 
		val(meta)
	
	output:
		path '*_general_stats.csv',          emit: stats
		path "${meta.id}_epicore.csv",       emit: final_epicore_csv
		path "length_distributions.html",    emit: length_dist
		path "epitope_intensity_hist.html",  emit: intensity_hist
		path "versions.yml",                 emit: versions   

	script:
		def args = task.ext.args ?: ''
		def prefix = task.ext.prefix ?: "${meta.id}"
		"""#!/bin/bash
		create_params_yaml.py \\
			--prefix $prefix \\
			$args

		epicore --reference_proteome $fasta --params_file ${prefix}_epicore_params.yml generate-epicore-csv --evidence_file $quantification_tsv

		(echo '<!DOCTYPE html><html><head><meta charset="UTF-8"><title>SVG</title></head><body>'; cat length_distributions.svg; echo '</body></html>') > length_distributions.html
		
		(echo '<!DOCTYPE html><html><head><meta charset="UTF-8"><title>SVG</title></head><body>'; cat epitope_intensity_hist.svg; echo '</body></html>') > epitope_intensity_hist.html

		cp pep_cores_mapping.csv ${prefix}_epicore.csv

		wc -l < epitopes.csv | awk '{print \$1 - 1}' > epicores.txt
		
		awk 'NR==1 {print \$0 ",# Epicores"; next} NR==2 {getline extra < "epicores.txt"; print \$0 "," extra}' $general_stats > modified_$general_stats
		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		    epicore: \$(epicore --version | grep 'epicore' | cut -d ' ' -f3)
		END_VERSIONS
		""" 

	stub:
		def args = task.ext.args ?: ''
		def prefix = task.ext.prefix ?: "${meta.id}"
		"""
		touch ${prefix}_general_stats.csv
		touch ${meta.id}_epicore.csv	
		touch length_distributions.html
		touch epitope_intensity_hist.html
		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		    epicore: \$(epicore --version | grep 'epicore' | cut -d ' ' -f3)
		END_VERSIONS		
		"""

}