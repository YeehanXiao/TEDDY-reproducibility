"""
Quantify read-level support for simulated intra-exonic TE-host breakpoints.

Evidence types:
1. Continuous CIGAR-defined alignment across a TE-host boundary.
2. Continuous full-span alignment across a short embedded TE segment.
3. Paired-anchor support linking TE-side and host-side regions.

Soft clipping is recorded only as an auxiliary diagnostic.
"""

import argparse
import csv
import glob
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import pysam


def get_nh(read):
  try:
  return int(read.get_tag("NH"))
except KeyError:
  return 1


def init_counter():
  return {
    "total": set(),
    "NH1": set(),
    "NHgt1": set()
  }


def add_read(counter, read):
  qname = read.query_name
nh = get_nh(read)

counter["total"].add(qname)

if nh == 1:
  counter["NH1"].add(qname)
elif nh > 1:
  counter["NHgt1"].add(qname)


def intersect_counter(a, b):
  return {
    "total": a["total"].intersection(b["total"]),
    "NH1": a["NH1"].intersection(b["NH1"]),
    "NHgt1": a["NHgt1"].intersection(b["NHgt1"])
  }


def parse_region(region_str):
  regions = []

if not region_str or region_str == "NA":
  return regions

for item in region_str.split(";"):
  item = item.strip()

if not item or item == "NA":
  continue

chrom, coords = item.split("|")
start, end = coords.split("-")

regions.append({
  "chrom": chrom,
  "start": int(start) - 1,
  "end": int(end)
})

return regions


def get_anchor_read_counter(
  bam,
  regions,
  min_overlap=15
):
  counter = init_counter()

for region in regions:
  try:
  for read in bam.fetch(
    region["chrom"],
    max(0, region["start"]),
    region["end"]
  ):
  if read.is_unmapped:
  continue

if read.get_overlap(
  region["start"],
  region["end"]
) >= min_overlap:
  add_read(counter, read)

except ValueError:
  pass

return counter


def check_continuous_boundary(
  blocks,
  pos_1based,
  anchor=10
):
  pos0 = pos_1based - 1

for block_start, block_end in blocks:
  left_ok = block_start <= (
    pos0 - anchor + 1
  )
right_ok = block_end >= (
  pos0 + anchor + 1
)

if left_ok and right_ok:
  return True

return False


def check_continuous_full_span(
  blocks,
  start_1based,
  end_1based,
  anchor=10
):
  start0 = start_1based - 1
end0 = end_1based - 1

for block_start, block_end in blocks:
  left_ok = block_start <= (
    start0 - anchor + 1
  )
right_ok = block_end >= (
  end0 + anchor + 1
)

if left_ok and right_ok:
  return True

return False


def has_softclip_near_boundary(
  read,
  pos_1based,
  tolerance=5
):
  if read.cigartuples is None:
  return False

has_softclip = any(
  operation == 4
  for operation, length in read.cigartuples
)

if not has_softclip:
  return False

ref_start_1based = read.reference_start + 1
ref_end_1based = read.reference_end

return (
  abs(ref_start_1based - pos_1based) <= tolerance
  or abs(ref_end_1based - pos_1based) <= tolerance
)


def get_cigar_read_counter(
  bam,
  chrom,
  start,
  end,
  evidence_type,
  anchor=10
):
  counter = init_counter()
softclip_counter = init_counter()

fetch_start = min(start, end) - anchor - 50
fetch_end = max(start, end) + anchor + 50

try:
  for read in bam.fetch(
    chrom,
    max(0, fetch_start - 1),
    fetch_end
  ):
  if read.is_unmapped:
  continue

blocks = read.get_blocks()

if not blocks:
  continue

is_supported = False

if evidence_type == "long_boundary_left":
  is_supported = check_continuous_boundary(
    blocks,
    start,
    anchor
  )

elif evidence_type == "long_boundary_right":
  is_supported = check_continuous_boundary(
    blocks,
    end,
    anchor
  )

elif evidence_type == "short_full_span":
  is_supported = check_continuous_full_span(
    blocks,
    start,
    end,
    anchor
  )

if is_supported:
  add_read(counter, read)

boundary_pos = (
  end
  if evidence_type == "long_boundary_right"
  else start
)

if has_softclip_near_boundary(
  read,
  boundary_pos
):
  add_read(softclip_counter, read)

except ValueError:
  pass

return counter, softclip_counter


def process_bam(
  bam_path,
  anchors,
  anchor_len=10,
  min_overlap=15
):
  bam = pysam.AlignmentFile(bam_path, "rb")

sample_name = os.path.basename(
  bam_path
).replace(".bam", "")

results = []

for row in anchors:
  te_regions = parse_region(
    row.get("TE_region", "")
  )
host_regions = parse_region(
  row.get("Host_region", "")
)

te_reads = get_anchor_read_counter(
  bam,
  te_regions,
  min_overlap
)

host_reads = get_anchor_read_counter(
  bam,
  host_regions,
  min_overlap
)

anchor_support = intersect_counter(
  te_reads,
  host_reads
)

chrom = row["chrom"]
start = int(row["anchor_start"])
end = int(row["anchor_end"])
evidence_type = row["evidence_type"]

cigar_support, softclip_diag = (
  get_cigar_read_counter(
    bam=bam,
    chrom=chrom,
    start=start,
    end=end,
    evidence_type=evidence_type,
    anchor=anchor_len
  )
)

results.append({
  "depth": sample_name,
  "bp_id": row["bp_id"],
  "transcript_id": row.get(
    "transcript_id",
    "NA"
  ),
  "gene_id": row.get(
    "gene_id",
    "NA"
  ),
  "gene_name": row.get(
    "gene_name",
    "NA"
  ),
  "TE_name": row.get(
    "TE_name",
    "NA"
  ),
  "evidence_type": evidence_type,
  
  "Anchor_total_reads": len(
    anchor_support["total"]
  ),
  "Anchor_NH1_reads": len(
    anchor_support["NH1"]
  ),
  "Anchor_NHgt1_reads": len(
    anchor_support["NHgt1"]
  ),
  "Anchor_supported": int(
    len(anchor_support["total"]) > 0
  ),
  "Anchor_unique_supported": int(
    len(anchor_support["NH1"]) > 0
  ),
  
  "CIGAR_total_reads": len(
    cigar_support["total"]
  ),
  "CIGAR_NH1_reads": len(
    cigar_support["NH1"]
  ),
  "CIGAR_NHgt1_reads": len(
    cigar_support["NHgt1"]
  ),
  "CIGAR_supported": int(
    len(cigar_support["total"]) > 0
  ),
  "CIGAR_unique_supported": int(
    len(cigar_support["NH1"]) > 0
  ),
  
  "softclip_near_boundary_reads": len(
    softclip_diag["total"]
  ),
  "softclip_near_boundary_NH1_reads": len(
    softclip_diag["NH1"]
  )
})

bam.close()
return results


def main():
  parser = argparse.ArgumentParser()

parser.add_argument(
  "-t",
  "--tsv",
  required=True,
  help="Breakpoint-anchor TSV"
)

parser.add_argument(
  "-d",
  "--bam_dir",
  required=True,
  help="Directory containing depth-specific BAM files"
)

parser.add_argument(
  "-o",
  "--out",
  required=True,
  help="Output TSV"
)

parser.add_argument(
  "--anchor_len",
  type=int,
  default=10
)

parser.add_argument(
  "--min_overlap",
  type=int,
  default=15
)

parser.add_argument(
  "--threads",
  type=int,
  default=8
)

args = parser.parse_args()

with open(args.tsv, "r") as handle:
  anchors = list(
    csv.DictReader(
      handle,
      delimiter="\t"
    )
  )

bam_files = sorted(
  glob.glob(
    os.path.join(
      args.bam_dir,
      "*.bam"
    )
  )
)

if not bam_files:
  raise FileNotFoundError(
    f"No BAM files found in {args.bam_dir}"
  )

print(
  f"Running breakpoint support analysis "
  f"on {len(anchors)} loci..."
)
print(f"BAM files: {len(bam_files)}")

all_rows = []

with ProcessPoolExecutor(
  max_workers=args.threads
) as executor:
  
  futures = {
    executor.submit(
      process_bam,
      bam,
      anchors,
      args.anchor_len,
      args.min_overlap
    ): bam
    for bam in bam_files
  }

for future in as_completed(futures):
  bam_path = futures[future]
all_rows.extend(future.result())
print(
  f"Done: {os.path.basename(bam_path)}"
)

out_dir = os.path.dirname(args.out)

if out_dir:
  os.makedirs(out_dir, exist_ok=True)

fields = [
  "depth",
  "bp_id",
  "transcript_id",
  "gene_id",
  "gene_name",
  "TE_name",
  "evidence_type",
  
  "Anchor_total_reads",
  "Anchor_NH1_reads",
  "Anchor_NHgt1_reads",
  "Anchor_supported",
  "Anchor_unique_supported",
  
  "CIGAR_total_reads",
  "CIGAR_NH1_reads",
  "CIGAR_NHgt1_reads",
  "CIGAR_supported",
  "CIGAR_unique_supported",
  
  "softclip_near_boundary_reads",
  "softclip_near_boundary_NH1_reads"
]

with open(
  args.out,
  "w",
  newline=""
) as handle:
  writer = csv.DictWriter(
    handle,
    fieldnames=fields,
    delimiter="\t"
  )
writer.writeheader()
writer.writerows(all_rows)

print(
  f"Breakpoint support table saved to "
  f"{args.out}"
)


if __name__ == "__main__":
  main()