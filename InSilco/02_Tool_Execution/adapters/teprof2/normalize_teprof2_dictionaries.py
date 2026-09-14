#!/usr/bin/env python2
from __future__ import print_function

import cPickle as pickle
import os
import shutil
import sys


def normalize(path):
    strand = "minus" if "minus" in os.path.basename(path) else "plus"

    with open(path, "rb") as handle:
        dic = pickle.load(handle)

    n_intron = 0
    missing_before = 0
    missing_exon_number = 0

    for chrom, transcripts in dic.items():
        for transcript, elements in transcripts.items():

            intron_keys = []

            for key, value in elements.items():
                fields = value.split(",")

                if len(fields) <= 8:
                    raise RuntimeError(
                        "Malformed dictionary record: {} {} {}".format(
                            chrom, transcript, value
                        )
                    )

                if fields[2] == "intron":
                    n_intron += 1
                    intron_keys.append(key)

                    if "intron_number " not in fields[8]:
                        missing_before += 1

                elif fields[2] == "exon":
                    if "exon_number " not in fields[8]:
                        missing_exon_number += 1

            # TEProf2 creates introns in increasing genomic order for both
            # strands, then reverses the numbering for minus-strand records
            # during downstream annotation.
            intron_keys = sorted(
                intron_keys,
                key=lambda x: (int(x[0]), int(x[1]))
            )

            for number, key in enumerate(intron_keys, 1):
                fields = elements[key].split(",")

                parts = [
                    part.strip()
                    for part in fields[8].split(";")
                    if part.strip()
                    and not part.strip().startswith("exon_number ")
                    and not part.strip().startswith("intron_number ")
                ]

                parts.append("intron_number {}".format(number))
                fields[8] = "; ".join(parts) + ";"

                elements[key] = ",".join(fields)

    if missing_exon_number:
        raise RuntimeError(
            "{} exon records lack exon_number in {}".format(
                missing_exon_number, path
            )
        )

    missing_after = 0

    for chrom, transcripts in dic.items():
        for transcript, elements in transcripts.items():
            for value in elements.values():
                fields = value.split(",")

                if (
                    len(fields) > 8
                    and fields[2] == "intron"
                    and "intron_number " not in fields[8]
                ):
                    missing_after += 1

    if missing_after:
        raise RuntimeError(
            "{} malformed intron records remain in {}".format(
                missing_after, path
            )
        )

    backup = path + ".before_intron_normalization"

    if not os.path.exists(backup):
        shutil.copy2(path, backup)

    tmp = path + ".tmp"

    with open(tmp, "wb") as handle:
        pickle.dump(dic, handle, protocol=2)

    os.rename(tmp, path)

    print(path)
    print("strand =", strand)
    print("intron records =", n_intron)
    print("missing intron_number before =", missing_before)
    print("missing intron_number after =", missing_after)


for dictionary in sys.argv[1:]:
    normalize(dictionary)
