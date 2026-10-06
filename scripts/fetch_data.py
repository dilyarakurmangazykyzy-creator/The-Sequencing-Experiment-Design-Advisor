#!/usr/bin/env python3
"""Retrieve public archived data, preserving bounded-prefix provenance."""
import argparse, csv, gzip, hashlib, io, json, pathlib, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIELDS = 'run_accession,study_accession,sample_accession,scientific_name,instrument_platform,instrument_model,library_layout,library_strategy,read_count,base_count,fastq_ftp,fastq_bytes,first_public'

def request(url):
    for attempt in range(4):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent':'SequencingAdvisor/1.0 academic reproducible retrieval'}), timeout=180)
        except Exception:
            if attempt == 3: raise
            time.sleep(2 ** attempt)

def report(accession):
    url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urllib.parse.urlencode(dict(accession=accession,result='read_run',fields=FIELDS,format='tsv'))
    with request(url) as response: text=response.read().decode()
    rows=list(csv.DictReader(io.StringIO(text),delimiter='\t'))
    if len(rows)!=1: raise ValueError(f'{accession}: expected exactly one resolved record, got {len(rows)}')
    return rows[0], url

def subset(url, target, count):
    start=time.perf_counter(); bases=0; got=0
    tmp=target.with_suffix(target.suffix+'.part')
    with request(url) as response, gzip.GzipFile(fileobj=response) as stream, tmp.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb',mtime=0,filename='') as out:
            for i in range(count):
                lines=[stream.readline() for _ in range(4)]
                if not lines[0]: break
                if not all(lines) or not lines[0].startswith(b'@') or not lines[2].startswith(b'+') or len(lines[1].strip())!=len(lines[3].strip()):
                    raise ValueError(f'invalid FASTQ record {i+1} at {url}')
                out.write(b''.join(lines)); bases+=len(lines[1].strip()); got+=1
    tmp.replace(target)
    return dict(file=target.name,source_url=url,records=got,bases=bases,selection='first N archived records; prefix acquisition is not a random whole-run sample',sha256=hashlib.sha256(target.read_bytes()).hexdigest(),elapsed_seconds=round(time.perf_counter()-start,3),archive_checksum_verified=False)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--short-pairs',type=int,default=400000);ap.add_argument('--long-reads',type=int,default=60000);ap.add_argument('--metadata-limit',type=int,default=100000);ap.add_argument('--out',type=pathlib.Path,default=ROOT/'data/raw');args=ap.parse_args()
    if min(args.short_pairs,args.long_reads,args.metadata_limit)<1: ap.error('sizes must be positive')
    args.out.mkdir(parents=True,exist_ok=True)
    rows=[]; sources=[]
    for acc in ['SRR1030394','SRR33602302']:
        row,url=report(acc);rows.append(row);sources.append(url);print('Resolved',acc,row['instrument_model'],flush=True)
    with (args.out/'selected_runs.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS.split(','),delimiter='\t');w.writeheader();w.writerows(rows)
    jobs=[]
    for i,url in enumerate(rows[0]['fastq_ftp'].split(';')):
        if url.endswith('_1.fastq.gz'): jobs.append(('https://'+url,args.out/'reads_illumina_R1.fastq.gz',args.short_pairs))
        if url.endswith('_2.fastq.gz'): jobs.append(('https://'+url,args.out/'reads_illumina_R2.fastq.gz',args.short_pairs))
    longs=rows[1]['fastq_ftp'].split(';')
    if len(longs)!=1: raise ValueError('long-read run must have one unambiguous FASTQ')
    jobs.append(('https://'+longs[0],args.out/'reads_nanopore.fastq.gz',args.long_reads))
    with ThreadPoolExecutor(max_workers=3) as pool:
        files=list(pool.map(lambda j:subset(*j),jobs))
    ref_url='https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3?download=false'
    with request(ref_url) as r: ref=r.read()
    if not ref.startswith(b'>') or len(ref)<4000000: raise ValueError('reference failed validation')
    # Preserve sequence accession in FASTA header. U00096.3 is GenBank counterpart of NC_000913.3.
    (args.out/'reference.fasta').write_bytes(ref)
    files.append(dict(file='reference.fasta',source_url=ref_url,accession='U00096.3',refseq_counterpart='NC_000913.3',sha256=hashlib.sha256(ref).hexdigest()))
    meta_fields='run_accession,study_accession,sample_accession,scientific_name,instrument_platform,instrument_model,library_layout,library_strategy,read_count,base_count,first_public'
    query='first_public>=2015-01-01 AND first_public<=2025-12-31 AND (instrument_platform="ILLUMINA" OR instrument_platform="OXFORD_NANOPORE" OR instrument_platform="PACBIO_SMRT")'
    meta_url='https://www.ebi.ac.uk/ena/portal/api/search?'+urllib.parse.urlencode(dict(result='read_run',query=query,fields=meta_fields,format='tsv',limit=args.metadata_limit))
    start=time.perf_counter(); meta_path=args.out/'metadata.tsv'
    with request(meta_url) as r,meta_path.open('wb') as f:
        while block:=r.read(1024*1024): f.write(block)
    with meta_path.open(encoding='utf-8') as f: n=sum(1 for _ in f)-1
    manifest=dict(retrieved_utc=datetime.now(timezone.utc).isoformat(),runs=rows,metadata_urls=sources,files=files,metadata=dict(source_url=meta_url,query=query,requested_limit=args.metadata_limit,rows=n,sha256=hashlib.sha256(meta_path.read_bytes()).hexdigest(),elapsed_seconds=round(time.perf_counter()-start,3),sampling='bounded unsorted ENA query; descriptive archive subset, not representative prevalence'),comparability='Same named MG1655 strain, different BioSamples/studies/years and WGS versus WGA. Observed difference cannot be attributed to platform alone.')
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__': main()
