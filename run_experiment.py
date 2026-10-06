#!/usr/bin/env python3
"""One-command reproducible real-read QC, minimap2 alignment and downsampling."""
import argparse, hashlib, json, pathlib, shutil, subprocess, sys, urllib.request

ROOT=pathlib.Path(__file__).resolve().parent

def run(args):
    print('Running:', ' '.join(map(str,args)),flush=True)
    subprocess.run(list(map(str,args)),cwd=ROOT,check=True)

def fetch_tools():
    directory=ROOT/'.cache/biowasm';directory.mkdir(parents=True,exist_ok=True)
    hashes=json.loads((ROOT/'config/biowasm_hashes.json').read_text())
    for extension,digest in hashes.items():
        local=directory/('minimap2.'+('cjs' if extension=='js' else extension))
        if not local.exists() or hashlib.sha256(local.read_bytes()).hexdigest()!=digest:
            url=f'https://biowasm.com/cdn/v3/minimap2/2.22/minimap2-simd.{extension}'
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 SequencingAdvisor/1.0'}),timeout=180) as response: content=response.read()
            if hashlib.sha256(content).hexdigest()!=digest: raise ValueError('Pinned tool checksum differs: '+url)
            local.write_bytes(content)
    return directory

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--pilot',action='store_true',help='reproduce full measured bounded pilot; fetches ~350 MB of compressed reads if absent');ap.add_argument('--node',default=shutil.which('node'));ap.add_argument('--force-fetch',action='store_true');args=ap.parse_args()
    if not args.node: ap.error('Install Node.js 18 or later, or supply --node PATH. The local advisor UI needs only Python.')
    data=ROOT/'data/raw' if args.pilot else ROOT/'data/example'
    out=ROOT/'results' if args.pilot else ROOT/'results/quick'
    if args.pilot and (args.force_fetch or not (data/'manifest.json').exists()): run([sys.executable,'scripts/fetch_data.py'])
    if not (data/'reference.fasta').exists(): ap.error('Included example data are missing; download the complete project bundle.')
    tools=fetch_tools();out.mkdir(parents=True,exist_ok=True);aln=out/'alignments';aln.mkdir(exist_ok=True)
    for platform,preset in [('illumina','sr'),('nanopore','map-ont')]:
        files=[data/'reference.fasta']+([data/'reads_illumina_R1.fastq.gz',data/'reads_illumina_R2.fastq.gz'] if platform=='illumina' else [data/'reads_nanopore.fastq.gz'])
        run([args.node,'--max-old-space-size=4096','scripts/align_wasm.cjs','--tools',tools,aln/(platform+'.sam'),'-ax',preset,'-t','1','--secondary=no',*files])
    run([sys.executable,'pipeline.py','pilot','--data',data,'--out',out,'--sam-illumina',aln/'illumina.sam','--sam-nanopore',aln/'nanopore.sam','--sam-aligner','minimap2_2.22_biowasm','--max-illumina','0','--max-nanopore','0'])
    if args.pilot:
        meta=data/'metadata_million.tsv'
        if not meta.exists():
            run([sys.executable,'scripts/metadata_mine.py','--fetch','--input',meta,'--out',out,'--query','first_public>=2015-01-01 AND first_public<=2025-12-31','--limit','1000000'])
        run([sys.executable,'scripts/metadata_mine.py','--input',meta,'--out',out])
    print('Completed. Measured outputs:',out,flush=True)

if __name__=='__main__':main()
