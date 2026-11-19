#!/bin/bash
#SBATCH --job-name=motifbreakR_SnpEff
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/motifbreakR_SnpEff_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/motifbreakR_SnpEff_%j.err
#SBATCH --partition=day
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Annotate ONLY HIGH_MEDIUM_SNPs in the existing gene-focused table.
# Usage: sbatch run_motifbreakR_SnpEff.sh
# Preflight: bash run_motifbreakR_SnpEff.sh --check
# Overrides: --root DIR --table TSV --vcf VCF.gz --reference FASTA --snpeff-db NAME
# Keep this and the evidence/plotting jobs sequential if they write the same table.
set -euo pipefail
CONFIG="${MOTIF_SNPEFF_CONFIG:-$HOME/donglab/pipelines/modules/motifbreakR_SnpEff/runtime.env}"
if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
elif [[ -n "${MOTIF_SNPEFF_CONFIG:-}" ]]; then
    echo "ERROR: configuration not found: $CONFIG" >&2; exit 1
fi
export PYTHONNOUSERSITE=1
"${ANNOTATION_PYTHON:-python3}" -u - "$@" <<'PYTHON'
import argparse, collections, csv, datetime, fcntl, gzip, hashlib, json
import math, os, re, shutil, subprocess, sys, tempfile
from pathlib import Path
csv.field_size_limit(100000000)

R_SOURCE = r'''
suppressPackageStartupMessages({
  library(motifbreakR)
  library(MotifDb)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BiocParallel)
})
args <- commandArgs(TRUE)
genome <- BSgenome.Hsapiens.UCSC.hg38
sources <- c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C")
motifs <- subset(MotifDb, dataSource %in% sources & organism == "Hsapiens")
if (!length(motifs)) stop("No human HOCOMOCO v11 A/B/C motifs in installed MotifDb")
cat("Human motifs:", length(motifs), "\n")
print(sessionInfo())
if (args[1] == "check") quit(status=0)
bed <- args[2]; outdir <- args[3]; workers <- as.integer(args[4]); threshold <- as.numeric(args[5])
input <- read.delim(bed, header=FALSE, stringsAsFactors=FALSE)
if (!nrow(input)) stop("Empty motif input")
if (any(!input$V1 %in% seqnames(genome))) stop("Input contig absent from hg38 BSgenome")
# Confirm the same reference alleles in the Bioconductor genome as in the FASTA.
expected <- vapply(strsplit(input$V4, ":", fixed=TRUE), function(z) z[3], character(1))
observed <- as.character(Biostrings::getSeq(genome, names=input$V1,
                                          start=input$V3, end=input$V3))
if (any(toupper(observed) != expected)) stop("REF mismatch against hg38 BSgenome")
variants <- snps.from.file(file=bed, search.genome=genome, format="bed",
                          check.unnamed.for.rsid=FALSE)
if (length(variants) != nrow(input)) stop("motifbreakR imported a different number of variants")
if (!setequal(as.character(variants$SNP_id), input$V4)) stop("motifbreakR changed variant identifiers")
backend <- if (workers > 1L) BiocParallel::MulticoreParam(workers=workers) else BiocParallel::SerialParam()
set.seed(1)
results <- motifbreakR(snpList=variants, pwmList=motifs, filterp=TRUE,
                      threshold=threshold, method="ic",
                      bkg=c(A=0.25,C=0.25,G=0.25,T=0.25), BPPARAM=backend)
saveRDS(results, file.path(outdir,"motifbreakR.results.rds"))
cols <- c("variant_key","chrom","start","end","strand","TF","motif_ID","source",
          "pctRef","pctAlt","scoreRef","scoreAlt","alleleDiff","effect")
if (length(results)) {
  m <- S4Vectors::mcols(results)
  take <- function(n) if (n %in% names(m)) as.character(m[[n]]) else rep(NA_character_,length(results))
  detail <- data.frame(variant_key=take("SNP_id"),chrom=as.character(seqnames(results)),
      start=start(results),end=end(results),strand=as.character(strand(results)),
      TF=take("geneSymbol"),motif_ID=take("providerId"),source=take("dataSource"),
      pctRef=take("pctRef"),pctAlt=take("pctAlt"),scoreRef=take("scoreRef"),
      scoreAlt=take("scoreAlt"),alleleDiff=take("alleleDiff"),effect=take("effect"),
      stringsAsFactors=FALSE)
  if (any(!detail$variant_key %in% input$V4)) stop("Unknown variant in motif results")
} else {
  detail <- as.data.frame(setNames(replicate(length(cols),character(),simplify=FALSE),cols))
}
write.table(detail,file.path(outdir,"motifbreakR.details.tsv"),sep="\t",quote=TRUE,
            row.names=FALSE,na="NA")
writeLines("completed",file.path(outdir,"motifbreakR.completed"))
'''

NEW_COLUMNS = ["SnpEff_by_SNP", "SnpEff_matching_gene_by_SNP", "motifbreakR_by_SNP",
               "functional_annotation_details"]

def fail(message):
    raise RuntimeError(message)

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def text_open(path):
    return gzip.open(path, 'rt', newline='') if str(path).endswith('.gz') else Path(path).open(newline='')

def read_table(path):
    with text_open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            fail(f'Empty or duplicate table headers: {path}')
        return reader.fieldnames, list(reader)

def write_table(path, fields, rows):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def selected_snps(rows):
    per_row = []
    for row in rows:
        values = sorted({x.strip() for x in row['HIGH_MEDIUM_SNPs'].split(',') if x.strip()})
        if not values:
            fail(f'Empty HIGH_MEDIUM_SNPs for gene {row.get("gene", "?")}')
        if any(not re.fullmatch(r'rs\d+', x) for x in values):
            fail('Expected comma-separated rsIDs in HIGH_MEDIUM_SNPs; no coordinate/allele guessing')
        per_row.append(values)
    return per_row, sorted({x for values in per_row for x in values})

def chrom_key(value):
    value = re.sub(r'^chr', '', str(value), flags=re.I).upper()
    return 'MT' if value == 'M' else value

def ucsc_chrom(value):
    return 'chrM' if chrom_key(value) == 'MT' else 'chr'+chrom_key(value)

def position(value):
    # Accept integer strings and exact integer decimal representations, never round.
    try:
        v = float(value)
        if not math.isfinite(v) or v < 1 or v != int(v): raise ValueError()
        return int(v)
    except (TypeError, ValueError, OverflowError):
        fail(f'Invalid GRCh38 position: {value!r}')

def lookup_positions(evidence, wanted):
    found = {}
    with text_open(evidence) as f:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'SNP', 'GWAS_chr', 'GWAS_pos'}
        if not required.issubset(reader.fieldnames or []): fail(f'Missing evidence columns: {required}')
        for row in reader:
            snp = row['SNP']
            if snp not in wanted: continue
            if snp in found: fail(f'Duplicate SNP in evidence table: {snp}')
            c = chrom_key(row['GWAS_chr'])
            if c not in {str(i) for i in range(1,23)} | {'X','Y','MT'}:
                fail(f'Unsupported chromosome for {snp}: {c}')
            found[snp] = dict(chrom=c, pos=position(row['GWAS_pos']),
                             effect_allele=row.get('GWAS_effect_allele','').upper(),
                             other_allele=row.get('GWAS_other_allele','').upper())
    missing = set(wanted)-set(found)
    if missing: fail('SNPs absent from evidence table: '+', '.join(sorted(missing)))
    return found

def plink_alleles(raw, locations):
    targets = {(v['chrom'],v['pos']) for v in locations.values()}
    found = collections.defaultdict(set)
    with Path(raw).open() as f:
        for token in f.readline().split()[6:]:
            m = re.fullmatch(r'(?:chr)?([^:]+):(\d+):([^:]+):([^_]+)_(.+)', token, flags=re.I)
            if not m: continue
            chrom,pos,ref,alt,counted = m.groups()
            key = (chrom_key(chrom),int(pos))
            if key in targets: found[key].add((ref.upper(),alt.upper(),counted.upper()))
    return found

def expected_variants(locations, alleles):
    result = {}; audit = []
    for snp,loc in locations.items():
        aa = alleles.get((loc['chrom'],loc['pos']),set())
        row = dict(SNP=snp,chrom=loc['chrom'],pos=loc['pos'])
        if len(aa) != 1:
            row['status'] = 'missing_or_ambiguous_PLINK_alleles'
        else:
            ref,alt,counted = next(iter(aa))
            row.update(REF=ref,ALT=alt)
            if ref not in 'ACGT' or alt not in 'ACGT' or len(ref)!=1 or len(alt)!=1 or ref==alt or counted not in (ref,alt):
                row['status'] = 'not_a_supported_biallelic_SNV'
            elif all(a in {'A','C','G','T'} for a in (loc['effect_allele'],loc['other_allele'])) and {loc['effect_allele'],loc['other_allele']} != {ref,alt}:
                row['status'] = 'GWAS_PLINK_allele_mismatch'
            else:
                key = f"{ucsc_chrom(loc['chrom'])}:{loc['pos']}:{ref}:{alt}"
                result[snp] = dict(loc, REF=ref,ALT=alt,variant_key=key)
                row['status'] = 'ready'
        audit.append(row)
    return result,audit

def fasta_index(path):
    result = {}
    with Path(str(path)+'.fai').open() as f:
        for line in f:
            cells = line.rstrip().split('\t')
            key = chrom_key(cells[0])
            if key in result: fail('Ambiguous FASTA contig aliases: '+key)
            result[key] = (cells[0], *map(int,cells[1:5]))
    return result

def reference_base(fasta, fai, chrom, pos):
    if chrom not in fai: fail('Contig absent from reference: '+chrom)
    name,length,offset,bases,width = fai[chrom]
    if not 1 <= pos <= length: fail(f'Position outside {name}: {pos}')
    zero = pos-1
    with Path(fasta).open('rb') as f:
        f.seek(offset+(zero//bases)*width+(zero%bases))
        return f.read(1).decode('ascii').upper()

def run(cmd, log, stdout=None, cwd=None):
    print('Running:', ' '.join(map(str,cmd)), flush=True)
    with Path(log).open('a') as err:
        if stdout:
            with Path(stdout).open('w') as out:
                subprocess.run(list(map(str,cmd)),stdout=out,stderr=err,cwd=cwd,check=True)
        else:
            subprocess.run(list(map(str,cmd)),stdout=err,stderr=subprocess.STDOUT,cwd=cwd,check=True)

def extract_candidates(args, expected, fai, folder):
    header = subprocess.check_output([args.bcftools,'view','-h',str(args.vcf)],text=True)
    contigs = {}
    for line in header.splitlines():
        m = re.match(r'##contig=<ID=([^,>]+)',line)
        if m:
            k = chrom_key(m[1])
            if k in contigs: fail('Ambiguous VCF contig aliases: '+k)
            contigs[k] = m[1]
    targets = sorted({(v['chrom'],v['pos']) for v in expected.values()})
    absent = {c for c,p in targets}-set(contigs)
    if absent: fail('Target chromosomes absent from VCF header: '+', '.join(sorted(absent)))
    regionfile = folder/'candidate_regions.tsv'
    with regionfile.open('w') as f:
        for c,p in targets: f.write(f'{contigs[c]}\t{p}\t{p}\n')
    raw = folder/'candidate_position_records.vcf'
    # Indexed seek to requested regions; drop all sample genotypes before writing.
    run([args.bcftools,'view','-G','-R',regionfile,'-Ov',args.vcf],folder/'extraction.log',stdout=raw)
    wanted = collections.defaultdict(set)
    for snp,v in expected.items(): wanted[(v['chrom'],v['pos'],v['REF'],v['ALT'])].add(snp)
    observed = set(); source_ids = collections.defaultdict(set)
    with raw.open() as f:
        for line in f:
            if line.startswith('#'): continue
            cells = line.rstrip('\n').split('\t')
            if len(cells)<8: fail('Malformed VCF record')
            c,p,ref = chrom_key(cells[0]),int(cells[1]),cells[3].upper()
            for alt in cells[4].upper().split(','):
                key = (c,p,ref,alt)
                if key in wanted:
                    observed.add(key);source_ids[key].add(cells[2])
    audit = []
    for snp,v in expected.items():
        key = (v['chrom'],v['pos'],v['REF'],v['ALT'])
        base = reference_base(args.reference,fai,v['chrom'],v['pos'])
        status = 'validated' if key in observed and base==v['REF'] else 'VCF_variant_missing' if key not in observed else 'FASTA_REF_mismatch'
        audit.append(dict(SNP=snp,variant_key=v['variant_key'],VCF_IDs=';'.join(sorted(source_ids[key])),status=status))
    write_table(folder/'variant_validation.tsv',['SNP','variant_key','VCF_IDs','status'],audit)
    if any(r['status']!='validated' for r in audit):
        fail('Some requested alleles could not be validated; see variant_validation.tsv. Table was not changed.')
    unique = {v['variant_key']:v for v in expected.values()}
    vcf = folder/'selected_snps.vcf';bed = folder/'selected_snps.hg38.bed'
    with vcf.open('w') as f,bed.open('w') as b:
        f.write('##fileformat=VCFv4.2\n##reference=GRCh38\n')
        for c in sorted({v['chrom'] for v in unique.values()}):
            f.write(f'##contig=<ID={ucsc_chrom(c)},length={fai[c][1]}>\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for key,v in sorted(unique.items()):
            c = ucsc_chrom(v['chrom']);p = v['pos']
            f.write(f"{c}\t{p}\t{key}\t{v['REF']}\t{v['ALT']}\t.\t.\t.\n")
            b.write(f'{c}\t{p-1}\t{p}\t{key}\t0\t+\n')
    print(f'Subset: {len(expected)} selected rsIDs, {len(unique)} unique REF/ALT variants. No whole-VCF annotation.',flush=True)
    return vcf,bed,unique

def parse_snpeff(path, valid_keys):
    fields = ['variant_key','allele','consequence','impact','gene','gene_id','feature_type',
              'transcript','biotype','rank','HGVS_c','HGVS_p','warnings','ANN_raw']
    details = []; seen = set()
    with Path(path).open() as f:
        for line in f:
            if line.startswith('#'): continue
            cells = line.rstrip('\n').split('\t');key=cells[2]
            if key not in valid_keys: fail('Unexpected variant in SnpEff output: '+key)
            if key in seen: fail('Duplicate variant in SnpEff output: '+key)
            seen.add(key)
            ann = next((x[4:] for x in cells[7].split(';') if x.startswith('ANN=')), '')
            for item in ann.split(',') if ann else []:
                a = item.split('|')
                if len(a)<16: fail('Unexpected SnpEff ANN format: '+item)
                if a[0] != cells[4]: fail('SnpEff annotation allele mismatch: '+item)
                details.append(dict(zip(fields,[key,*a[:10],a[10],a[15],item])))
    if seen != set(valid_keys): fail('SnpEff output is missing selected variants')
    return fields,details

def as_json(value):
    return json.dumps(value,ensure_ascii=False,separators=(',',':'))

def merge_annotations(rows, per_row, expected, snpeff, motifs, folder):
    se = collections.defaultdict(list);mb = collections.defaultdict(list)
    for r in snpeff:
        se[r['variant_key']].append({k:r[k] for k in ['consequence','impact','gene','gene_id','transcript','HGVS_c','HGVS_p','warnings']})
    for r in motifs:
        if r['variant_key'] not in {v['variant_key'] for v in expected.values()}:
            fail('Unexpected motif variant: '+r['variant_key'])
        mb[r['variant_key']].append({k:r.get(k,'') for k in ['TF','motif_ID','effect','pctRef','pctAlt','scoreRef','scoreAlt','alleleDiff','chrom','start','end','strand','source']})
    answer = []
    for row,ids in zip(rows,per_row):
        new = dict(row);snp_data={};gene_data={};motif_data={}
        gene = row.get('gene','').upper()
        for snp in ids:
            key=expected[snp]['variant_key']
            effects=se[key]
            snp_data[snp]={'variant':key,'status':'annotated' if effects else 'no_ANN_returned','effects':effects}
            matches=[r for r in effects if r['gene'].upper()==gene or re.sub(r'\.\d+$','',r['gene_id']).upper()==re.sub(r'\.\d+$','',gene)]
            gene_data[snp]={'status':'matched' if matches else 'no_matching_gene_annotation','effects':matches}
            motif_data[snp]={'variant':key,'status':'motif_hits' if mb[key] else 'no_hits_at_selected_threshold',
                             'hits':mb[key]}
        new.update(SnpEff_by_SNP=as_json(snp_data),SnpEff_matching_gene_by_SNP=as_json(gene_data),
                   motifbreakR_by_SNP=as_json(motif_data),functional_annotation_details=str(folder))
        answer.append(new)
    return answer

def atomic_merge(table, fields, before, after, digest, folder):
    if sha256(table) != digest: fail('Input table changed while annotation ran; refusing to overwrite it')
    if len(before)!=len(after): fail('Row count changed during merge')
    original_fields = [c for c in fields if c not in NEW_COLUMNS]
    for old,new in zip(before,after):
        if any(old[c]!=new[c] for c in original_fields): fail('An existing field changed during annotation')
    columns=original_fields+NEW_COLUMNS
    backup=folder/(table.name+'.before_annotation')
    shutil.copy2(table,backup)
    fd,temp=tempfile.mkstemp(prefix='.'+table.name+'.',suffix='.tmp',dir=table.parent)
    os.close(fd)
    try:
        write_table(temp,columns,after)
        f2,r2=read_table(temp)
        if f2!=columns or r2!=after: fail('Saved table verification failed')
        os.chmod(temp,table.stat().st_mode & 0o777)
        with open(temp,'rb') as f: os.fsync(f.fileno())
        if sha256(table)!=digest: fail('Table changed before final replacement')
        os.replace(temp,table)
    finally:
        if Path(temp).exists(): Path(temp).unlink()
    return backup

def arguments():
    p=argparse.ArgumentParser(description='Annotate only HIGH_MEDIUM_SNPs and merge into the existing gene-focused table.')
    root=Path(os.environ.get('TARGET_ALS_ROOT',str(Path.home()/'donglab/data/target_ALS')))
    p.add_argument('--root',type=Path,default=root)
    p.add_argument('--table',type=Path)
    p.add_argument('--evidence',type=Path)
    p.add_argument('--vcf',type=Path)
    p.add_argument('--reference',type=Path,default=Path(os.environ.get('REFERENCE_FASTA',str(Path.home()/'donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'))))
    p.add_argument('--plink-raw',type=Path)
    p.add_argument('--snpeff-db',default=os.environ.get('SNPEFF_DB',''))
    modulehome=Path.home()/'donglab/pipelines/modules'
    p.add_argument('--snpeff-jar',type=Path,default=Path(os.environ.get('SNPEFF_JAR',str(modulehome/'snpEff/snpEff.jar'))))
    p.add_argument('--snpeff-config',type=Path,default=Path(os.environ.get('SNPEFF_CONFIG',str(modulehome/'snpEff/snpEff.config'))))
    p.add_argument('--snpeff-data',type=Path,default=Path(os.environ.get('SNPEFF_DATA',str(modulehome/'snpEff/data'))))
    p.add_argument('--bcftools',default=os.environ.get('BCFTOOLS_BIN','bcftools'))
    p.add_argument('--java',default=os.environ.get('JAVA_BIN','java'))
    p.add_argument('--rscript',default=os.environ.get('RSCRIPT_BIN','Rscript'))
    p.add_argument('--threads',type=int,default=int(os.environ.get('SLURM_CPUS_PER_TASK','4')))
    p.add_argument('--motif-threshold',type=float,default=1e-4)
    p.add_argument('--check',action='store_true',help='Check selected SNPs, installed tools/references and R dependencies; no annotation or table edits.')
    a=p.parse_args()
    a.root=a.root.expanduser().resolve()
    out=a.root/'QTL/ALS_variant_validation_evidence'
    a.table=(a.table or out/'gene_focused_interesting_SNPs_from_QTL.tsv').expanduser().resolve()
    a.evidence=(a.evidence or out/'ALS_GWAS_variant_validation_evidence.ALL.tsv.gz').expanduser().resolve()
    a.vcf=(a.vcf or Path(os.environ.get('TARGET_ALS_VCF',str(a.root/'QTL/target_ALS_imputed_filtered_joint.vcf.gz')))).expanduser().resolve()
    a.plink_raw=(a.plink_raw or a.root/'QTL/plink/joint_all_chrs_matrixEQTL.raw').expanduser().resolve()
    for key in ['reference','snpeff_jar','snpeff_config','snpeff_data']:
        setattr(a,key,getattr(a,key).expanduser().resolve())
    if a.threads<1 or a.threads>int(os.environ.get('SLURM_CPUS_PER_TASK',str(a.threads))):
        p.error('--threads must be positive and no greater than allocated Slurm CPUs')
    if not 0<a.motif_threshold<1: p.error('--motif-threshold must be between 0 and 1')
    return a

def main():
    a=arguments()
    for path in [a.table,a.evidence,a.vcf,a.reference,a.plink_raw]:
        if not path.is_file(): fail(f'Missing input: {path}')
    fields,rows=read_table(a.table)
    if not {'gene','HIGH_MEDIUM_SNPs'}.issubset(fields): fail('Table requires gene and HIGH_MEDIUM_SNPs columns')
    per_row,ids=selected_snps(rows)
    if not ids:
        print('No HIGH/MEDIUM SNPs; nothing to annotate or merge.');return
    print(f'{len(rows)} gene rows; {len(ids)} unique requested rsIDs. Input VCF: {a.vcf}',flush=True)
    missing=[]
    for command in [a.bcftools,a.java,a.rscript]:
        if not shutil.which(command): missing.append('Executable not available: '+command)
    for path in [a.snpeff_jar,a.snpeff_config,Path(str(a.reference)+'.fai')]:
        if not path.is_file(): missing.append('Missing installed file: '+str(path))
    if not any(Path(str(a.vcf)+suffix).is_file() for suffix in ['.tbi','.csi']):
        missing.append('VCF requires an existing .tbi or .csi index')
    if not a.snpeff_db:
        missing.append('Set SNPEFF_DB or --snpeff-db to the installed GRCh38 database name')
    elif 'GRCh38' not in a.snpeff_db and 'hg38' not in a.snpeff_db:
        missing.append('SnpEff database name must identify GRCh38/hg38')
    elif not (a.snpeff_data/a.snpeff_db/'snpEffectPredictor.bin').is_file():
        missing.append('SnpEff database not installed: '+str(a.snpeff_data/a.snpeff_db))
    if missing: fail('Preflight needs setup:\n  '+'\n  '.join(missing))
    if a.check:
        with tempfile.TemporaryDirectory(prefix='motif_snpeff_check_') as tmp:
            rfile=Path(tmp)/'motifbreakR.R';rfile.write_text(R_SOURCE)
            run([a.rscript,'--vanilla',rfile,'check'],Path(tmp)/'R.log')
            print((Path(tmp)/'R.log').read_text(),flush=True)
        locations=lookup_positions(a.evidence,set(ids));expected,audit=expected_variants(locations,plink_alleles(a.plink_raw,locations))
        if len(expected)!=len(ids): fail('PLINK/GWAS allele check failed: '+as_json([r for r in audit if r['status']!='ready']))
        fai=fasta_index(a.reference)
        for snp,v in expected.items():
            if reference_base(a.reference,fai,v['chrom'],v['pos'])!=v['REF']: fail('FASTA REF mismatch: '+snp)
        print('Preflight passed; no predictions run and no tables changed.');return
    lockpath=Path(str(a.table)+'.motif_SnpEff.lock')
    with lockpath.open('a') as lock:
        try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: fail('Another annotation job is already updating this table')
        digest=sha256(a.table)
        # Re-read under the lock and bind results to exactly this table snapshot.
        fields,rows=read_table(a.table);per_row,ids=selected_snps(rows)
        stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        base=a.table.parent/'motifbreakR_SnpEff';base.mkdir(parents=True,exist_ok=True)
        folder=Path(tempfile.mkdtemp(prefix=stamp+'_job'+os.environ.get('SLURM_JOB_ID','local')+'_',dir=base))
        print('Annotation details:',folder,flush=True)
        write_table(folder/'requested_rsIDs.tsv',['SNP'],[dict(SNP=x) for x in ids])
        locations=lookup_positions(a.evidence,set(ids))
        expected,audit=expected_variants(locations,plink_alleles(a.plink_raw,locations))
        write_table(folder/'allele_mapping_audit.tsv',['SNP','chrom','pos','REF','ALT','status'],audit)
        if len(expected)!=len(ids): fail('Unresolved/unsupported variants; see allele_mapping_audit.tsv. Input table unchanged.')
        fai=fasta_index(a.reference)
        vcf,bed,unique=extract_candidates(a,expected,fai,folder)
        rfile=folder/'run_motifbreakR.R';rfile.write_text(R_SOURCE)
        run([a.rscript,'--vanilla',rfile,'check'],folder/'motifbreakR.log')
        # Annotation is offline; no database auto-download or usage reporting.
        snpeff_vcf=folder/'selected_snps.SnpEff.vcf'
        command=[a.java,'-Xmx12g','-jar',a.snpeff_jar,'ann','-noLog','-nodownload','-noStats',
                 '-c',a.snpeff_config,'-dataDir',a.snpeff_data,a.snpeff_db,vcf]
        run(command,folder/'SnpEff.log',stdout=snpeff_vcf,cwd=folder)
        se_fields,se=parse_snpeff(snpeff_vcf,set(unique))
        write_table(folder/'SnpEff.details.tsv',se_fields,se)
        run([a.rscript,'--vanilla',rfile,'run',bed,folder,str(a.threads),str(a.motif_threshold)],folder/'motifbreakR.log')
        if not (folder/'motifbreakR.completed').is_file(): fail('motifbreakR did not complete')
        _,motifs=read_table(folder/'motifbreakR.details.tsv')
        after=merge_annotations(rows,per_row,expected,se,motifs,folder)
        summary=[]
        for snp,v in sorted(expected.items()):
            hits=[r for r in motifs if r['variant_key']==v['variant_key']]
            summary.append(dict(SNP=snp,variant_key=v['variant_key'],
                SnpEff_annotations=sum(r['variant_key']==v['variant_key'] for r in se),
                motif_hits=len(hits),strong_motif_hits=sum(r['effect']=='strong' for r in hits)))
        write_table(folder/'per_SNP_annotation_summary.tsv',list(summary[0]),summary)
        metadata=dict(start_UTC=stamp,arguments={k:str(v) for k,v in vars(a).items()},
            input_table_sha256=digest,evidence_file=dict(path=str(a.evidence),size=a.evidence.stat().st_size,mtime=a.evidence.stat().st_mtime),
            input_vcf=dict(path=str(a.vcf),size=a.vcf.stat().st_size,mtime=a.vcf.stat().st_mtime),
            selected_rsIDs=len(ids),unique_allelic_variants=len(unique),motifs='Human HOCOMOCOv11 core A/B/C',
            motif_method='ic',motif_threshold=a.motif_threshold,
            interpretation='Sequence-motif predictions only, not measured TF occupancy or disease causality. No scores or lead SNPs changed.')
        for name,cmd in [('bcftools',[a.bcftools,'--version']),('java',[a.java,'-version']),
                         ('SnpEff',[a.java,'-jar',str(a.snpeff_jar),'-version'])]:
            v=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            metadata[name+'_version']=v.stdout.strip()
        (folder/'run_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
        backup=atomic_merge(a.table,fields,rows,after,digest,folder)
        print(f'DONE: {a.table}\nBackup: {backup}\nDetails: {folder}',flush=True)
        print('Only four annotation columns added/refreshed. Scores, priority, lead SNPs and existing rows preserved.',flush=True)

if __name__=='__main__':
    try: main()
    except (RuntimeError,OSError,ValueError,subprocess.CalledProcessError) as error:
        print('ERROR:',error,file=sys.stderr,flush=True)
        sys.exit(1)
PYTHON
