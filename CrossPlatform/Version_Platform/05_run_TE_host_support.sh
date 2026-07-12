#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash 03_run_TE_host_support.sh \
    <TE_host_pair_meta.tsv> <TE_host_eval_anchors.bed> <out_dir> \
    <threads_per_job> <max_parallel_jobs> \
    <sample1=bam1> [sample2=bam2 ...]

Example:
  bash 03_run_TE_host_support.sh \
    output/TE_host_pair_meta.tsv \
    output/TE_host_eval_anchors.bed \
    output/nanopore_support 4 3 \
    control_1=input/control_1.bam \
    control_2=input/control_2.bam \
    control_3=input/control_3.bam
USAGE
}

if [[ $# -lt 6 ]]; then
  usage
  exit 1
fi

meta="$1"
eval_bed="$2"
out_dir="$3"
threads="$4"
max_jobs="$5"
shift 5
sample_specs=("$@")

for cmd in samtools bedtools awk xargs; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Missing command: $cmd" >&2; exit 1; }
done

[[ -s "$meta" ]] || { echo "Missing meta file: $meta" >&2; exit 1; }
[[ -s "$eval_bed" ]] || { echo "Missing anchor BED: $eval_bed" >&2; exit 1; }

mkdir -p \
  "$out_dir/bam_subset" \
  "$out_dir/read_bed" \
  "$out_dir/intersect" \
  "$out_dir/per_read" \
  "$out_dir/summary" \
  "$out_dir/logs"

col_idx() {
  local file="$1"
  local name="$2"
  awk -v target="$name" 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i==target){print i; exit}}' "$file"
}

C_PAIR=$(col_idx "$meta" pair_id)
C_TX=$(col_idx "$meta" transcript_id)
C_GENE=$(col_idx "$meta" gene_name)
C_EVAL=$(col_idx "$meta" eval_type)
C_HOSTTYPE=$(col_idx "$meta" host_anchor_type)
C_TELEN=$(col_idx "$meta" te_len)
C_HOSTLEN=$(col_idx "$meta" host_len)
C_PAIRLEN=$(col_idx "$meta" pair_len)

for x in C_PAIR C_TX C_GENE C_EVAL C_HOSTTYPE C_TELEN C_HOSTLEN C_PAIRLEN; do
  [[ -n "${!x}" ]] || { echo "Required column missing from meta: $x" >&2; exit 1; }
done

jobs_file="$out_dir/jobs.tsv"
: > "$jobs_file"
for spec in "${sample_specs[@]}"; do
  sample="${spec%%=*}"
  bam="${spec#*=}"
  [[ "$sample" != "$bam" ]] || { echo "Invalid sample specification: $spec" >&2; exit 1; }
  [[ -s "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }
  printf '%s\t%s\n' "$sample" "$bam" >> "$jobs_file"
done

run_intersect_one() {
  local sample="$1"
  local bam="$2"

  echo "[$(date)] START alignment extraction: $sample"

  samtools view -@ "$threads" -b -L "$eval_bed" "$bam" \
    > "$out_dir/bam_subset/${sample}.TE_host_anchors.bam" \
    2> "$out_dir/logs/${sample}.samtools_subset.log"

  bedtools bamtobed -bed12 \
    -i "$out_dir/bam_subset/${sample}.TE_host_anchors.bam" \
    > "$out_dir/read_bed/${sample}.reads.bed12" \
    2> "$out_dir/logs/${sample}.bamtobed.log"

  bedtools intersect -split \
    -a "$out_dir/read_bed/${sample}.reads.bed12" \
    -b "$eval_bed" \
    -wo \
    > "$out_dir/intersect/${sample}.read_vs_TE_host_anchors.wo.tsv"

  echo "[$(date)] DONE alignment extraction: $sample"
}

export -f run_intersect_one
export threads eval_bed out_dir
xargs -P "$max_jobs" -n 2 bash -c 'run_intersect_one "$1" "$2"' _ < "$jobs_file"

make_per_read_one() {
  local sample="$1"
  local intersect_file="$out_dir/intersect/${sample}.read_vs_TE_host_anchors.wo.tsv"
  local per_read_file="$out_dir/per_read/${sample}.per_read_TE_host_anchor_coverage.tsv"

  [[ -s "$intersect_file" ]] || { echo "Missing intersect file: $intersect_file" >&2; exit 1; }
  echo "[$(date)] START per-read summary: $sample"

  awk \
    -v sample="$sample" \
    -v c_pair="$C_PAIR" -v c_tx="$C_TX" -v c_gene="$C_GENE" \
    -v c_eval="$C_EVAL" -v c_hosttype="$C_HOSTTYPE" \
    -v c_telen="$C_TELEN" -v c_hostlen="$C_HOSTLEN" -v c_pairlen="$C_PAIRLEN" '
    BEGIN{FS=OFS="\t"}
    NR==FNR{
      if(FNR==1){next}
      pair=$c_pair
      pair_tx[pair]=$c_tx
      pair_gene[pair]=$c_gene
      pair_eval[pair]=$c_eval
      pair_hosttype[pair]=$c_hosttype
      te_len[pair]=$c_telen+0
      host_len[pair]=$c_hostlen+0
      pair_len[pair]=$c_pairlen+0
      next
    }
    {
      read_id=$4
      split($16,a,"|")
      pair=a[1]
      role=a[4]
      host_type=a[5]
      ov=$19+0
      key=pair SUBSEP read_id

      pair_seen[key]=pair
      read_seen[key]=read_id
      host_type_seen[key]=host_type
      if(role=="TE_anchor") te_cov[key]+=ov
      if(role=="host_anchor") host_cov[key]+=ov
    }
    END{
      print "sample","pair_id","transcript_id","gene_name","eval_type","host_anchor_type","read_id","TE_cov_bases","host_cov_bases","TE_len","host_len","pair_len","TE_cov_frac","host_cov_frac","pair_cov_frac","TE_host_junction_10bp","TE_host_anchor75"
      for(key in pair_seen){
        pair=pair_seen[key]
        if(!(pair in te_len)) continue
        if(te_len[pair]<=0 || host_len[pair]<=0) continue

        te=te_cov[key]+0
        host=host_cov[key]+0
        if(te>te_len[pair]) te=te_len[pair]
        if(host>host_len[pair]) host=host_len[pair]

        te_frac=te/te_len[pair]
        host_frac=host/host_len[pair]
        pair_frac=(te+host)/(te_len[pair]+host_len[pair])
        junction=(te>=10 && host>=10 ? "Yes" : "-")
        anchor75=(te_frac>=0.75 && host_frac>=0.75 ? "Yes" : "-")

        print sample,pair,pair_tx[pair],pair_gene[pair],pair_eval[pair],host_type_seen[key],read_seen[key],te,host,te_len[pair],host_len[pair],pair_len[pair],te_frac,host_frac,pair_frac,junction,anchor75
      }
    }
  ' "$meta" "$intersect_file" > "$per_read_file"

  echo "[$(date)] DONE per-read summary: $sample"
}

export -f make_per_read_one
export meta out_dir C_PAIR C_TX C_GENE C_EVAL C_HOSTTYPE C_TELEN C_HOSTLEN C_PAIRLEN
cut -f1 "$jobs_file" | xargs -P "$max_jobs" -n 1 bash -c 'make_per_read_one "$1"' _

mapfile -t samples < <(cut -f1 "$jobs_file")
merged_per_read="$out_dir/per_read/combined.per_read_TE_host_anchor_coverage.tsv"
cat "$out_dir/per_read/${samples[0]}.per_read_TE_host_anchor_coverage.tsv" > "$merged_per_read"
for sample in "${samples[@]:1}"; do
  tail -n +2 "$out_dir/per_read/${sample}.per_read_TE_host_anchor_coverage.tsv" >> "$merged_per_read"
done

pair_out="$out_dir/summary/pair_TE_host_read_support.tsv"
tx_out="$out_dir/summary/tx_TE_host_read_support.tsv"
rate_out="$out_dir/summary/TE_host_read_support_rates.tsv"

awk \
  -v c_pair="$C_PAIR" -v c_tx="$C_TX" -v c_gene="$C_GENE" \
  -v c_eval="$C_EVAL" -v c_hosttype="$C_HOSTTYPE" \
  -v c_telen="$C_TELEN" -v c_hostlen="$C_HOSTLEN" -v c_pairlen="$C_PAIRLEN" '
  BEGIN{FS=OFS="\t"}
  NR==FNR{
    if(FNR==1){next}
    pair=$c_pair
    pair_tx[pair]=$c_tx
    pair_gene[pair]=$c_gene
    pair_eval[pair]=$c_eval
    pair_hosttype[pair]=$c_hosttype
    te_len[pair]=$c_telen+0
    host_len[pair]=$c_hostlen+0
    pair_len[pair]=$c_pairlen+0
    all_pair[pair]=1
    next
  }
  FNR==1{next}
  {
    pair=$2
    all_pair[pair]=1
    n_reads_any[pair]++
    if(($13+0)>max_te[pair]) max_te[pair]=$13+0
    if(($14+0)>max_host[pair]) max_host[pair]=$14+0
    if(($15+0)>max_pair[pair]) max_pair[pair]=$15+0
    if($16=="Yes") n_junction[pair]++
    if($17=="Yes") n_anchor75[pair]++
  }
  END{
    print "pair_id","transcript_id","gene_name","eval_type","host_anchor_type","TE_len","host_len","pair_len","n_reads_any","max_TE_cov_frac","max_host_cov_frac","max_pair_cov_frac","n_TE_host_junction_reads_10bp","TE_host_junction_support_10bp","n_TE_host_anchor75_reads","TE_host_anchor75_support"
    for(pair in all_pair){
      print pair,pair_tx[pair],pair_gene[pair],pair_eval[pair],pair_hosttype[pair],te_len[pair],host_len[pair],pair_len[pair],n_reads_any[pair]+0,max_te[pair]+0,max_host[pair]+0,max_pair[pair]+0,n_junction[pair]+0,((n_junction[pair]+0)>0?"Yes":"-"),n_anchor75[pair]+0,((n_anchor75[pair]+0)>0?"Yes":"-")
    }
  }
' "$meta" "$merged_per_read" > "$pair_out"

awk \
  -v c_tx="$C_TX" -v c_gene="$C_GENE" -v c_eval="$C_EVAL" '
  BEGIN{FS=OFS="\t"}
  NR==FNR{
    if(FNR==1){next}
    tx=$c_tx
    all_tx[tx]=1
    tx_gene[tx]=$c_gene
    tx_eval[tx]=$c_eval
    next
  }
  FNR==1{next}
  {
    tx=$2
    all_tx[tx]=1
    tx_gene[tx]=$3
    tx_eval[tx]=$4
    n_pairs[tx]++
    n_reads[tx]+=$9
    if(($10+0)>max_te[tx]) max_te[tx]=$10+0
    if(($11+0)>max_host[tx]) max_host[tx]=$11+0
    if(($12+0)>max_pair[tx]) max_pair[tx]=$12+0
    n_junction[tx]+=$13
    n_anchor75[tx]+=$15
    if($14=="Yes") junction_yes[tx]=1
    if($16=="Yes") anchor75_yes[tx]=1
  }
  END{
    print "transcript_id","gene_name","eval_type","n_pairs","n_reads_any","max_TE_cov_frac","max_host_cov_frac","max_pair_cov_frac","n_TE_host_junction_reads_10bp","TE_host_junction_support_10bp","n_TE_host_anchor75_reads","TE_host_anchor75_support"
    for(tx in all_tx){
      print tx,tx_gene[tx],tx_eval[tx],n_pairs[tx]+0,n_reads[tx]+0,max_te[tx]+0,max_host[tx]+0,max_pair[tx]+0,n_junction[tx]+0,((junction_yes[tx]+0)>0?"Yes":"-"),n_anchor75[tx]+0,((anchor75_yes[tx]+0)>0?"Yes":"-")
    }
  }
' "$meta" "$pair_out" > "$tx_out"

awk '
  BEGIN{FS=OFS="\t"}
  NR==1{next}
  {
    n_all++
    if($3=="multi_exon"){
      n_multi++
      if($10=="Yes") j_multi++
      if($12=="Yes") a_multi++
    }
    if($3=="single_exon"){
      n_single++
      if($10=="Yes") j_single++
      if($12=="Yes") a_single++
    }
    if($10=="Yes") j_all++
    if($12=="Yes") a_all++
  }
  END{
    print "metric","numerator","denominator","rate_pct"
    print "TE_host_junction_support_10bp_multi_exon",j_multi+0,n_multi+0,100*(j_multi+0)/(n_multi+0)
    print "TE_host_junction_support_10bp_single_exon",j_single+0,n_single+0,100*(j_single+0)/(n_single+0)
    print "TE_host_junction_support_10bp_all",j_all+0,n_all+0,100*(j_all+0)/(n_all+0)
    print "TE_host_anchor75_support_multi_exon",a_multi+0,n_multi+0,100*(a_multi+0)/(n_multi+0)
    print "TE_host_anchor75_support_single_exon",a_single+0,n_single+0,100*(a_single+0)/(n_single+0)
    print "TE_host_anchor75_support_all",a_all+0,n_all+0,100*(a_all+0)/(n_all+0)
  }
' "$tx_out" > "$rate_out"

cat "$rate_out"
echo "Generated:"
echo "  $merged_per_read"
echo "  $pair_out"
echo "  $tx_out"
echo "  $rate_out"
