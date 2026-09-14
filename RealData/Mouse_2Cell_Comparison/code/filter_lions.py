#!/usr/bin/env python3

import argparse, re, subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig


COLS = [
    "transcriptID","exonRankInTranscript","repeatName","coordinates",
    "ER_Interaction","IsExonic","ExonsOverlappingWithRepeat",
    "ER","DR","DE","DD","Total","Chromosome",
    "EStart","EEnd","RStart","REnd","EStrand","RStrand",
    "RepeatRank","UpExonStart","UpExonEnd","UpThread","DownThread",
    "ExonInGene"
]

COV = [
    "ExonRPKM","ExonMax",
    "UpExonRPKM","UpExonMax",
    "RepeatRPKM","RepeatMaxCoverage",
    "UpstreamRepeatRPKM","UpstreamRepeatMaxCoverage"
]

# Real-data LIONS derives these coverage fields from the strand bigWigs used by
# its native workflow. The simulated-data adapter instead derives peaks from BAM.


def find_raw(repdir):
    x = sorted(repdir.rglob("*.lions.tsv"))
    if len(x) != 1:
        raise RuntimeError(f"{repdir}: expected 1 *.lions.tsv, found {len(x)}")
    return x[0]


def read_raw(path):
    with path.open() as handle:
        first = handle.readline().split("\t")[0]
    if first == "transcriptID":
        x = pd.read_csv(path, sep="\t", dtype=str)
        x.columns = COLS
    else:
        x = pd.read_csv(path, sep="\t", header=None, names=COLS, dtype=str)
    if x.shape[1] != 25:
        raise RuntimeError(f"{path}: {x.shape[1]} columns, expected 25")
    return x


def mapped_reads(bam, samtools):
    s = subprocess.check_output(
        [samtools, "flagstat", str(bam)],
        text=True
    )
    for line in s.splitlines():
        m = re.match(r"^(\d+)\s+\+\s+\d+\s+mapped\s+\(", line)
        if m:
            return int(m.group(1))
    raise RuntimeError(f"Cannot parse mapped reads: {bam}")


class BW:
    def __init__(self, fwd, rev):
        self.f = pyBigWig.open(str(fwd))
        self.r = pyBigWig.open(str(rev))
        self.fc = self.f.chroms()
        self.rc = self.r.chroms()
        self.cache = {}

    def close(self):
        self.f.close()
        self.r.close()

    def calc(self, chrom, start, end):
        try:
            start, end = int(float(start)), int(float(end))
        except Exception:
            return 0.0, 0.0

        key = (str(chrom), start, end)
        if key in self.cache:
            return self.cache[key]

        chrom = str(chrom)
        choices = [chrom]
        if chrom.startswith("chr"):
            choices.append(chrom[3:])
        else:
            choices.append("chr" + chrom)

        c = next(
            (z for z in choices if z in self.fc and z in self.rc),
            None
        )

        if c is None:
            self.cache[key] = (0.0, 0.0)
            return self.cache[key]

        end = min(end, self.fc[c], self.rc[c])
        start = max(0, start)

        if end <= start:
            self.cache[key] = (0.0, 0.0)
            return self.cache[key]

        a = np.asarray(self.f.values(c, start, end), dtype=float)
        b = np.asarray(self.r.values(c, start, end), dtype=float)

        v = (
            np.abs(np.nan_to_num(a, nan=0.0)) +
            np.abs(np.nan_to_num(b, nan=0.0))
        )

        out = (
            float(v.mean()) if len(v) else 0.0,
            float(v.max()) if len(v) else 0.0
        )
        self.cache[key] = out
        return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bench", required=True, type=Path)
    p.add_argument("--bw-dir", required=True, type=Path)
    p.add_argument("--bam-dir", required=True, type=Path)
    p.add_argument("--lions-dir", required=True, type=Path)
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--samtools", default="samtools")
    args = p.parse_args()

    work = args.bench / "LIONS" / "final_filter"
    work.mkdir(parents=True, exist_ok=True)

    chim = args.lions_dir / "scripts/ChimericReadTool/chimSort.R"
    qc = []

    for i in range(1, 5):
        sample = f"2cellrep{i}"
        rawfile = find_raw(args.bench / "LIONS" / sample)
        bam = args.bam_dir / f"{sample}.bam"
        fwd = args.bw_dir / f"{sample}_forward.bw"
        rev = args.bw_dir / f"{sample}_reverse.bw"

        print(f"\n===== {sample} =====", flush=True)
        print("raw:", rawfile, flush=True)

        x = read_raw(rawfile)

        for c in [
            "EStart","EEnd","RStart","REnd",
            "UpExonStart","UpExonEnd"
        ]:
            x[c] = pd.to_numeric(x[c], errors="coerce")

        bw = BW(fwd, rev)

        vals = {c: [] for c in COV}

        for row in x.itertuples(index=False):
            emean, emax = bw.calc(
                row.Chromosome, row.EStart, row.EEnd
            )
            umean, umax = bw.calc(
                row.Chromosome, row.UpExonStart, row.UpExonEnd
            )
            rmean, rmax = bw.calc(
                row.Chromosome, row.RStart, row.REnd
            )

            plus = str(row.EStrand) in ("1", "1.0", "+")
            if plus:
                us, ue = row.RStart - 50, row.RStart
            else:
                us, ue = row.REnd, row.REnd + 50

            trmean, trmax = bw.calc(row.Chromosome, us, ue)

            for c, v in zip(
                COV,
                [emean,emax,umean,umax,rmean,rmax,trmean,trmax]
            ):
                vals[c].append(v)

        bw.close()

        for c in COV:
            x[c] = vals[c]

        lcsv = work / f"{sample}.pc.lcsv"
        lion = work / f"{sample}.lion"

        x[COLS + COV].to_csv(
            lcsv, sep="\t", index=False
        )

        mapped = mapped_reads(bam, args.samtools)

        cmd = [
            args.rscript,
            str(chim),
            str(lcsv),
            str(lion),
            str(mapped),
            "3","10","10","1","0.1","2","1.5"
        ]

        subprocess.run(cmd, check=True)

        out = pd.read_csv(lion, sep="\t", dtype=str)
        final_genes = out["transcriptID"].str.split(":").str[0].nunique()
        effective_sc_reads = max(3, round(mapped / 20_000_000))
        qc.append({
            "replicate": f"rep{i}",
            "raw_interactions": len(x),
            "mapped_reads": mapped,
            "effective_scREADS": effective_sc_reads,
            "calls": len(out),
            "genes": final_genes,
        })
        print(
            "raw =", len(x),
            "mapped =", mapped,
            "scREADS =", effective_sc_reads,
            "final =", len(out),
            flush=True
        )

    pd.DataFrame(qc).to_csv(
        work / "lions_filter_qc.tsv",
        sep="\t",
        index=False,
    )

    print("\nLIONS_FILTER_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
