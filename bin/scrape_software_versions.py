#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/mhcquant': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'Comet': ['comet.params.new', r"comet_version (\S+)"],
    'Percolator': ['v_percolator.txt', r"version (\S+)"],
    'MHCFlurry': ['v_mhcflurry.txt', r"mhcflurry (\S+)"],
    'OpenMS': ['v_openms.txt', r"Version: (\S+)"],
}
results = OrderedDict()
results['nf-core/mhcquant'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Comet'] = '<span style="color:#999999;\">N/A</span>'
results['Percolator'] = '<span style="color:#999999;\">N/A</span>'
results['OpenMS'] = '<span style="color:#999999;\">N/A</span>'
results['MHCFlurry']= '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Remove software set to false in results
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/mhcquant Software Versions'
section_href: 'https://github.com/nf-core/mhcquant'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
