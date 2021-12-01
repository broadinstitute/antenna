version 1.0

task bam_to_fastq {
  input {
    File input_bam

    String docker = "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(size(input_bam, "Gi") * 10) + 10
    Int preemptible = 3
  }
  String picard_exec = "/usr/gitc/picard.jar"
  String output_fastq_1_filename = "output_1.fq"
  String output_fastq_2_filename = "output_2.fq"
  command <<<
    java -jar ~{picard_exec} SamToFastq I=~{input_bam} FASTQ=~{output_fastq_1_filename} SECOND_END_FASTQ=~{output_fastq_2_filename} 
  >>>
  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  output {
    File fastq1 = "~{output_fastq_1_filename}"
    File fastq2 = "~{output_fastq_2_filename}"
  }
}

task count_bam_reads {
  input {
    File input_bam
    String flags = ""

    String docker = "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(size(input_bam, "Gi") * 10) + 10
    Int preemptible = 3
  }

  String counts_filename = "counts.txt"

  command <<<
    samtools view ~{flags} -c ~{input_bam} > ~{counts_filename}
  >>>

  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    Int count = read_int(counts_filename)
  }
}

task bwa_align {
  input {
    File reference_bundle
    File fastq_r1
    File fastq_r2 

    String docker = "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438"
    Int machine_mem_mb = 32000
    Int cpu = 16
    Int disk = ceil(size(reference_bundle, "Gi") * 2 + size(fastq_r1, "Gi") * 2 + size(fastq_r2, "Gi") *2) + 10
    Int preemptible = 3
  }
  String output_aligned_bam_filename = "aligned.bam"
  command <<<
    set -euo pipefail

    tar --no-same-owner -xvf ~{reference_bundle}

    /usr/gitc/bwa mem -Y -t ~{cpu} ref/nCoV-2019.reference.fasta ~{fastq_r1} ~{fastq_r2} | samtools sort -@ 2 -m 8G - | samtools view -bh > ~{output_aligned_bam_filename}

    samtools index ~{output_aligned_bam_filename}
  >>>
  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks:  "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  output {
    File output_aligned_bam = "~{output_aligned_bam_filename}"
    File output_aligned_bam_index = "~{output_aligned_bam_filename}.bai"
  }
}

task antenna_tag {
  input {
      File inbam
      File inbam_index
      String outbam_name
      
      String docker = "quay.io/nbarkas_1/antenna:0.0.10"
      Int machine_mem_mb = 8192
      Int cpu = 1
      Int disk = ceil(size(inbam, "Gi") * 4) + 10
      Int preemptible = 3
  }

  String antenna_tag_exec = "antenna_tag_reads"
  
  command<<<
     ~{antenna_tag_exec} --bam ~{inbam} --outbam ~{outbam_name} 
  >>>
  
  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File outbam = "~{outbam_name}"
  }
}

task bam_sort_index {
  input {
    File bamfile

    String docker = "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438"
    Int machine_mem_mb = 32000
    Int cpu = 4
    Int disk = ceil(size(bamfile, "Gi") * 2 + 10)
    Int preemptible = 3
  }

  String output_aligned_bam_filename = "aligned.bam"

  command <<<
    set -euo pipefail

    samtools sort -@ 4 -m 6G ~{bamfile} > ~{output_aligned_bam_filename}

    samtools index ~{output_aligned_bam_filename}
  >>>

  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks:  "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File output_aligned_bam = "~{output_aligned_bam_filename}"
    File output_aligned_bam_index = "~{output_aligned_bam_filename}.bai"
  }
}

task antenna_count {
  input {
    File inbam
    File inbamindex
    File bedfile
    String outcsv_filename
    Int cutoff = "50"
    
    String docker = "quay.io/nbarkas_1/antenna:0.0.10"
    Int machine_mem_mb = 8192
    Int cpu = 1
    Int disk = ceil(size(inbam, "Gi") * 4) + 10
    Int preemptible = 3
  }
  
  String antenna_count_exec = "antenna_count_reads"

  command<<<
    ~{antenna_count_exec} --bam ~{inbam} --bed ~{bedfile} --outcsv ~{outcsv_filename} --cutoff ~{cutoff}
  >>>
  
  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File outcsv = "~{outcsv_filename}"
  } 
}

workflow antenna {
    input {
         File input_bam
         File ref_bundle
         File orf_locations
    }

    String version = "antenna_v0.0.10"

    call count_bam_reads as input_count {
      input:
        input_bam = input_bam
    }

    call bam_to_fastq {
      input:
        input_bam = input_bam
    }

    call bwa_align {
      input:
        reference_bundle = ref_bundle,
        fastq_r1 = bam_to_fastq.fastq1,
        fastq_r2 = bam_to_fastq.fastq2
    }

    call count_bam_reads as output_count {
      input:
        input_bam = bwa_align.output_aligned_bam,
        flags = " -F 260 "
    }
    
    call antenna_tag {
      input:
        inbam = bwa_align.output_aligned_bam,
        inbam_index = bwa_align.output_aligned_bam_index,
        outbam_name = "output.bam"
    }

    call bam_sort_index {
      input:
        bamfile = antenna_tag.outbam
    }
    
    call antenna_count {
      input:
         inbam = bam_sort_index.output_aligned_bam,
         inbamindex = bam_sort_index.output_aligned_bam_index,
         bedfile = orf_locations,
         outcsv_filename = "counts.csv"
    }

   output {
      Int input_read_count = input_count.count
      Int output_read_count = output_count.count
      File output_counts = antenna_count.outcsv
      File aligned_bam = bwa_align.output_aligned_bam
      String pipeline_version = "~{version}"
   }
}

