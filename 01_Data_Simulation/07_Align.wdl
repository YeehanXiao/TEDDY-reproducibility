version 1.0

workflow AlignSimulationReads {
  input {
    Array[String] samples
    String fastq_dir
    String star_index
    String bam_dir
    String rmdup_dir
    String sortbam_dir
    String picard_jar
    Int threads
  }

  scatter (sample in samples) {
    call STARbam {
      input:
        sampleName = sample,
        fastq_dir = fastq_dir,
        star_index = star_index,
        outdir = bam_dir,
        threads = threads
    }

    call MarkDuplicates {
      input:
        sampleName = sample,
        bampath = STARbam.bampath,
        rmdup_dir = rmdup_dir,
        picard_jar = picard_jar
    }

    call SortBam {
      input:
        sampleName = sample,
        bampath = MarkDuplicates.new_bampath,
        sortbam_dir = sortbam_dir,
        picard_jar = picard_jar
    }

    call BuildBamIndex {
      input:
        sampleName = sample,
        bampath = SortBam.sortbam_dir,
        picard_jar = picard_jar
    }
  }
}

task STARbam {
  input {
    String sampleName
    String fastq_dir
    String star_index
    String outdir
    Int threads
  }

  File read1 = "~{fastq_dir}/~{sampleName}.R1.fastq.gz"
  File read2 = "~{fastq_dir}/~{sampleName}.R2.fastq.gz"

  command <<<
    mkdir -p ~{outdir}

    STAR \
      --runThreadN ~{threads} \
      --genomeDir ~{star_index} \
      --readFilesCommand zcat \
      --readFilesIn ~{read1} ~{read2} \
      --outFileNamePrefix ~{outdir}/~{sampleName} \
      --outSAMtype BAM SortedByCoordinate
  >>>

  output {
    File bam = "~{outdir}/~{sampleName}Aligned.sortedByCoord.out.bam"
    String bampath = outdir
  }
}

task MarkDuplicates {
  input {
    String sampleName
    String bampath
    String rmdup_dir
    String picard_jar
  }

  command <<<
    mkdir -p ~{rmdup_dir}

    java -jar ~{picard_jar} MarkDuplicates \
      --REMOVE_DUPLICATES TRUE \
      -I ~{bampath}/~{sampleName}Aligned.sortedByCoord.out.bam \
      -O ~{rmdup_dir}/~{sampleName}.bam \
      -M ~{rmdup_dir}/~{sampleName}.metrics.txt \
      --ASSUME_SORT_ORDER coordinate
  >>>

  output {
    File rmdup_bam = "~{rmdup_dir}/~{sampleName}.bam"
    File metrics = "~{rmdup_dir}/~{sampleName}.metrics.txt"
    String new_bampath = rmdup_dir
  }
}

task SortBam {
  input {
    String sampleName
    String bampath
    String sortbam_dir
    String picard_jar
  }

  command <<<
    mkdir -p ~{sortbam_dir}

    java -jar ~{picard_jar} SortSam \
      -I ~{bampath}/~{sampleName}.bam \
      -O ~{sortbam_dir}/~{sampleName}.bam \
      --SORT_ORDER coordinate
  >>>

  output {
    File sorted_bam = "~{sortbam_dir}/~{sampleName}.bam"
    String sortbam_dir = sortbam_dir
  }
}

task BuildBamIndex {
  input {
    String sampleName
    String bampath
    String picard_jar
  }

  command <<<
    java -jar ~{picard_jar} BuildBamIndex \
      -I ~{bampath}/~{sampleName}.bam
  >>>

  output {
    File bai = "~{bampath}/~{sampleName}.bai"
  }
}