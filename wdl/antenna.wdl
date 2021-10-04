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

task bwa_align {
  input {
    File reference_bundle
    File fastq_r1
    File fastq_r2 

    String docker = "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438"
    Int machine_mem_mb = 32000
    Int cpu = 16
    Int disk = ceil(size(reference_bundle) * 10 + size(fastq_r1) * 2 + size(fastq_r2) *2) + 10
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

task antenna_task {
  input {
    File input_bam
    File input_bam_index
    File orf_locations 

    String docker = "quay.io/nbarkas_1/antenna:0.0.1"
    Int machine_mem_mb = 16
    Int cpu = 1
    Int disk = ceil(size(input_bam) * 10) + 10
    Int preemptible = 3
  }
  String antenna_exec = "/root/tools/antenna.py"
  String output_counts_filename = "counts.csv"
  command<<<
    ~{antenna_exec} --bam ~{input_bam} --output-counts ~{output_counts_filename} --orf-bed ~{orf_locations}
  >>>
  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  output {
    File counts = "~{output_counts_filename}"
  }
}    

workflow antenna {
    input {
         File input_bam
         File ref_bundle
         File orf_locations
    }

    String version = "antenna_v0.0.1"

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

    call antenna_task {
      input:
        input_bam = bwa_align.output_aligned_bam,
        input_bam_index = bwa_align.output_aligned_bam_index,
        orf_locations = orf_locations
    }

   output {
      File output_counts = antenna_task.counts
      String pipeline_version = "~{version}"
   }
}

