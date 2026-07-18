#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# Script: 02_quantify_chimeric_rpm.py
# Description: Quantifies transcript-level chimeric junctions by identifying reads 
#              where one end anchors within the TE (>= 15bp) and the mate maps 
#              outside the coordinates of that incorporated MTA sequence.
# ==============================================================================

import pysam
import glob
import os
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_single_bam(bam_path, bed_regions, min_overlap=15):
    """
    Process a single BAM file to count reads spanning out of TE regions.
    Returns the sample name, a dictionary of counts per region, and total mapped reads.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    
    # Record total mapped reads for downstream RPM normalization
    total_mapped = bam.mapped
    counts_dict = {}
    
    for region in bed_regions:
        chrom = region['chrom']
        start = region['start']
        end = region['end']
        name = region['name']
        
        valid_read_names = set()
        
        try:
            for read in bam.fetch(chrom, start, end):
                # Ensure the read is paired, the mate is mapped, and both map to the same chromosome
                if read.is_paired and not read.mate_is_unmapped and read.reference_id == read.next_reference_id:
                    
                    # Constraint 1: The read must overlap the TE region by at least the minimum required base pairs (e.g., 15bp)
                    if read.get_overlap(start, end) >= min_overlap:
                        
                        mate_start = read.next_reference_start
                        
                        # Constraint 2: The mate read must map outside the physical coordinates of the TE
                        if mate_start >= end or mate_start < start:
                            valid_read_names.add(read.query_name)
                            
            counts_dict[name] = len(valid_read_names)
            
        except ValueError:
            # Chromosome not found in BAM file
            counts_dict[name] = 0
            
    bam.close()
    return sample_name, counts_dict, total_mapped


def main():
    parser = argparse.ArgumentParser(description="Quantify transcript-level TE-chimeric junctions.")
    parser.add_argument("-b", "--bed", required=True, help="Input BED file containing TE overlap regions.")
    parser.add_argument("-d", "--bam_dir", required=True, help="Directory containing sorted BAM files.")
    parser.add_argument("-o", "--out_matrix", default="TE_spanning_counts_matrix.txt", help="Output count matrix file name.")
    parser.add_argument("-s", "--out_depth", default="sample_depths.txt", help="Output sequencing depth file name.")
    parser.add_argument("-m", "--min_overlap", type=int, default=15, help="Minimum overlap (bp) required within the TE (default: 15).")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads for parallel processing (default: 8).")
    args = parser.parse_args()

    # 1. Parse the input BED file
    bed_regions = []
    with open(args.bed, 'r') as bed:
        for line in bed:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            # Assumes standard BED format where column 4 is the locus/TE name
            bed_regions.append({
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'name': parts[3] 
            })

    # 2. Locate all BAM files
    bam_files = sorted(glob.glob(os.path.join(args.bam_dir, "*.bam")))
    if not bam_files:
        print(f"Error: No BAM files found in {args.bam_dir}")
        return

    all_results = {}
    sample_depths = {}

    print(f"Starting analysis on {len(bam_files)} BAM files...")
    print(f"Minimum TE overlap threshold: {args.min_overlap} bp")

    # 3. Execute parallel processing
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_single_bam, bam, bed_regions, args.min_overlap): bam for bam in bam_files}
        
        for future in as_completed(futures):
            try:
                sample_name, counts_dict, total_mapped = future.result()
                all_results[sample_name] = counts_dict
                sample_depths[sample_name] = total_mapped
                print(f"Processed sample: {sample_name:<15} | Total Mapped Reads: {total_mapped:,}")
            except Exception as e:
                print(f"Error processing {futures[future]}: {e}")

    # 4. Write Count Matrix
    print("\nWriting count matrix...")
    sample_names = sorted(list(all_results.keys()))
    
    with open(args.out_matrix, 'w') as out:
        header = ["Region_Name"] + sample_names
        out.write("\t".join(header) + "\n")
        
        for region in bed_regions:
            region_name = region['name']
            row = [region_name]
            for sample in sample_names:
                row.append(str(all_results[sample].get(region_name, 0)))
            out.write("\t".join(row) + "\n")

    # 5. Write Sequencing Depths
    with open(args.out_depth, 'w') as f_depth:
        f_depth.write("Sample\tTotal_Mapped_Reads\n")
        for sample in sample_names:
            f_depth.write(f"{sample}\t{sample_depths[sample]}\n")
            
    print(f"Analysis complete.")
    print(f"Count matrix saved to: {args.out_matrix}")
    print(f"Sequencing depths saved to: {args.out_depth}")

if __name__ == '__main__':
    main()