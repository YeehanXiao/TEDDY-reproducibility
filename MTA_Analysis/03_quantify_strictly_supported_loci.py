import pysam
import glob
import os
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed


def parse_region_string(region_string):
    """
    Parse collapsed neighbor regions.

    Expected format:
    chr1|4840956-4841132|upstream_exon|ENSMUST...|Gene|rank9

    Multiple regions are separated by ';'.

    Input coordinates are 1-based inclusive from R.
    Output coordinates are 0-based half-open for pysam.
    """
    regions = []

    if region_string is None:
        return regions

    region_string = str(region_string).strip()

    if region_string in ["", "NA", "<NA>", "nan", "NaN"]:
        return regions

    for item in region_string.split(";"):
        item = item.strip()
        if item == "":
            continue

        parts = item.split("|")
        if len(parts) < 3:
            continue

        chrom = parts[0]
        start_end = parts[1]
        region_type = parts[2]

        if "-" not in start_end:
            continue

        try:
            start_1based, end_1based = start_end.split("-", 1)
            start_1based = int(start_1based)
            end_1based = int(end_1based)
        except ValueError:
            continue

        if end_1based < start_1based:
            continue

        regions.append({
            "chrom": chrom,
            "start": start_1based - 1,
            "end": end_1based,
            "region_type": region_type,
            "raw": item
        })

    return regions


def read_locus_support_table(support_file):
    """
    Read locus-level support table.

    Required columns:
      chrom
      chromStart
      chromEnd
      name
      left_neighbor_regions
      right_neighbor_regions

    chrom/chromStart/chromEnd correspond to the MTA TE-overlap locus.
    Neighbor regions contain all possible transcript-context regions collapsed by ';'.
    """
    loci = []

    with open(support_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            te = {
                "chrom": row["chrom"],
                "start": int(row["chromStart"]),  # BED 0-based
                "end": int(row["chromEnd"]),
                "name": row["name"],
                "locus_id": row.get("locus_id", row["name"]),
                "TE_name": row.get("TE_name", "")
            }

            neighbors = []
            neighbors.extend(parse_region_string(row.get("left_neighbor_regions", "")))
            neighbors.extend(parse_region_string(row.get("right_neighbor_regions", "")))

            loci.append({
                "te": te,
                "neighbors": neighbors
            })

    return loci


def get_read_names_in_region(bam, region, min_overlap=15, min_mapq=0):
    """
    Return query names of primary reads overlapping region by >= min_overlap.
    CIGAR-aware through read.get_overlap().
    """
    read_names = set()

    try:
        for read in bam.fetch(region["chrom"], region["start"], region["end"]):
            if read.is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.get_overlap(region["start"], region["end"]) >= min_overlap:
                read_names.add(read.query_name)

    except ValueError:
        pass

    return read_names


def process_single_bam(
    bam_path,
    loci,
    min_te_overlap=15,
    min_neighbor_overlap=15,
    min_mapq=0
):
    bam = pysam.AlignmentFile(bam_path, "rb")
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    total_mapped = bam.mapped

    counts_dict = {}
    detail_rows = []

    for item in loci:
        te = item["te"]
        neighbors = item["neighbors"]
        locus_name = te["name"]

        if len(neighbors) == 0:
            counts_dict[locus_name] = 0
            continue

        te_read_names = get_read_names_in_region(
            bam=bam,
            region=te,
            min_overlap=min_te_overlap,
            min_mapq=min_mapq
        )

        if len(te_read_names) == 0:
            counts_dict[locus_name] = 0
            continue

        neighbor_read_names_all = set()
        neighbor_hit_map = {}

        for nb in neighbors:
            nb_read_names = get_read_names_in_region(
                bam=bam,
                region=nb,
                min_overlap=min_neighbor_overlap,
                min_mapq=min_mapq
            )

            if len(nb_read_names) == 0:
                continue

            neighbor_read_names_all.update(nb_read_names)

            for read_name in nb_read_names:
                if read_name not in neighbor_hit_map:
                    neighbor_hit_map[read_name] = []
                neighbor_hit_map[read_name].append(nb["raw"])

        supported_read_names = te_read_names.intersection(neighbor_read_names_all)

        counts_dict[locus_name] = len(supported_read_names)

        for read_name in supported_read_names:
            detail_rows.append({
                "sample": sample_name,
                "locus_id": te["locus_id"],
                "region_name": locus_name,
                "TE_name": te["TE_name"],
                "read_name": read_name,
                "neighbor_regions": ";".join(sorted(set(neighbor_hit_map.get(read_name, []))))
            })

    bam.close()

    return sample_name, counts_dict, total_mapped, detail_rows


if __name__ == "__main__":

    support_file = "/mnt/datadisk/xiaoyihan/manuals/junction/MTA_locus_neighbor_support_table.tsv"
    bam_dir = "/mnt/datadisk/xiaoyihan/TEchimeric/RNA/chenfei/sortBam"

    output_matrix = "MTA_locus_neighbor_spanning_counts_matrix_TE15_neighbor15.txt"
    depth_output = "sample_depths_TE15_neighbor15.txt"
    detail_output = "MTA_locus_neighbor_spanning_detail_TE15_neighbor15.tsv"

    min_te_overlap = 15
    min_neighbor_overlap = 15
    min_mapq = 0
    max_cores = 25

    loci = read_locus_support_table(support_file)
    bam_files = sorted(glob.glob(os.path.join(bam_dir, "*.bam")))

    print("Using query-name intersection logic")
    print(f"Support file: {support_file}")
    print(f"Loaded MTA loci: {len(loci)}")
    print(f"Min TE overlap: {min_te_overlap} bp")
    print(f"Min neighbor overlap: {min_neighbor_overlap} bp")
    print(f"Min MAPQ: {min_mapq}")
    print(f"Threads: {max_cores}")
    print(f"BAM files: {len(bam_files)}")
    print("Processing...\n")

    all_results = {}
    sample_depths = {}
    all_detail_rows = []

    with ProcessPoolExecutor(max_workers=max_cores) as executor:
        futures = {
            executor.submit(
                process_single_bam,
                bam,
                loci,
                min_te_overlap,
                min_neighbor_overlap,
                min_mapq
            ): bam
            for bam in bam_files
        }

        for future in as_completed(futures):
            bam_path = futures[future]

            try:
                sample_name, counts_dict, total_mapped, detail_rows = future.result()

                all_results[sample_name] = counts_dict
                sample_depths[sample_name] = total_mapped
                all_detail_rows.extend(detail_rows)

                supported_loci = sum(v > 0 for v in counts_dict.values())
                total_support = sum(counts_dict.values())

                print(
                    f"done {sample_name:<15} | mapped reads: {total_mapped:,} | "
                    f"supported loci: {supported_loci:,} | support reads: {total_support:,}"
                )

            except Exception as e:
                print(f"error processing {bam_path}: {e}")

    sample_names = sorted(all_results.keys())

    print("\nWriting count matrix...")

    with open(output_matrix, "w") as out:
        header = ["Region_Name"] + sample_names
        out.write("\t".join(header) + "\n")

        for item in loci:
            region_name = item["te"]["name"]
            row = [region_name]

            for sample in sample_names:
                row.append(str(all_results[sample].get(region_name, 0)))

            out.write("\t".join(row) + "\n")

    print("Writing depth file...")

    with open(depth_output, "w") as f_depth:
        f_depth.write("Sample\tTotal_Mapped_Reads\n")
        for sample in sample_names:
            f_depth.write(f"{sample}\t{sample_depths[sample]}\n")

    print("Writing read-level detail table...")

    detail_fields = [
        "sample",
        "locus_id",
        "region_name",
        "TE_name",
        "read_name",
        "neighbor_regions"
    ]

    with open(detail_output, "w") as f_detail:
        writer = csv.DictWriter(f_detail, fieldnames=detail_fields, delimiter="\t")
        writer.writeheader()

        for row in all_detail_rows:
            writer.writerow(row)

    print(f"\nCount matrix saved to: {output_matrix}")
    print(f"Depth file saved to: {depth_output}")
    print(f"Read-level detail saved to: {detail_output}")
