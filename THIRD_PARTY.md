# Third-party data and software

The real read fixture is a bounded excerpt of public INSDC/ENA runs SRR1030394 and SRR33602302. Reference accession U00096.3 is the GenBank counterpart of NC_000913.3. Source URLs and retrieval provenance are included in result manifests. The source runs are not generated or simulated by this project.

Alignment uses minimap2 2.22-r1101 (Heng Li, MIT license), compiled to WebAssembly by the Biowasm project (MIT license). Assets are retrieved on demand from the official Biowasm CDN and verified with SHA-256. This distribution does not bundle the binary assets. Cite Li (2018), doi:10.1093/bioinformatics/bty191, and Biowasm https://github.com/biowasm/biowasm .

Python/Node.js and any optional conda bioinformatics tools retain their upstream licenses. The full workflow lists its versions in `workflow/environment.yaml`.
