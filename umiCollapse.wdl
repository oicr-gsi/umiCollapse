version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem
import "imports/pull_bamQC.wdl" as bamQC

struct GenomeResources {
    String bwaMem_runBwaMem_bwaRef
    String bwaMem_runBwaMem_modules
}

struct FastqInputs {
    File fastq1
    File fastq2
    String readGroups
}

workflow umiCollapse {
    input {
        String umiList
        String outputPrefix = "output"
        Array[FastqInputs] fastqInputs
        String pattern1 
        String pattern2
        String reference
        Boolean doBamQC = false
        Boolean provisionBam = true
        String mode
    }

    Map[String, GenomeResources] resources = {
    "hg19": {
        "bwaMem_runBwaMem_bwaRef": "$HG19_BWA_INDEX_ROOT/hg19_random.fa",
        "bwaMem_runBwaMem_modules": "samtools/1.9 bwa/0.7.12 hg19-bwa-index/0.7.12"
    },
    "hg38": {
        "bwaMem_runBwaMem_bwaRef": "$HG38_BWA_INDEX_ROOT/hg38_random.fa",
        "bwaMem_runBwaMem_modules": "samtools/1.9 bwa/0.7.12 hg38-bwa-index/0.7.12"
    },
    "mm10": {
        "bwaMem_runBwaMem_bwaRef": "$MM10_BWA_INDEX_ROOT/mm10.fa", 
        "bwaMem_runBwaMem_modules": "samtools/1.9 bwa/0.7.12 mm10-bwa-index/0.7.12"
    }
    }

    parameter_meta {
        umiList: "Reference file with valid UMIs"
        outputPrefix: "Specifies the prefix of output files"
        fastqInputs: "Array of fastq structs containing reads and readgroups"
        pattern1: "UMI pattern 1"
        pattern2: "UMI pattern 2"
        reference: "Name and version of reference genome"
        doBamQC: "Enable/disable bamQC process"
        provisionBam: "Enable/disable provision out umi deduplicated bam file"
        mode: "running mode for the workflow, only allow value 'lane_level' and 'call_ready'"
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "The incorporation of Unique Molecular Indices (UMIs) into sequenced reads allows for more accurate identification of PCR duplicates. This workflow extracts UMIs from the reads in fastq files, into the sequence identifier line. UMI identification is based on a known location and pattern and with the ability to match to a list of expected sequences. Reads are then aligned to a reference genome with bwa, and then the aligned sequence file is collapsed to remove duplicates with um-Tools"
        dependencies: [
            {
                name: "barcodex-rs/0.1.2",
                url: "https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz"
            },
            {
                name: "rust/1.2",
                url: "https://www.rust-lang.org/tools/install"
            },
            {
                name: "umi-tools/1.1.1",
                url: "https://github.com/CGATOxford/UMI-tools/archive/1.1.1.tar.gz"
            },
            {
                name: "bwa/0.7.12",
                url: "https://github.com/lh3/bwa/archive/0.7.12.tar.gz"
            },
            {
                name: "samtools/1.9",
                url: "https://github.com/samtools/samtools/archive/0.1.19.tar.gz"
            },
            {
                name: "python/3.6",
                url: "https://www.python.org/downloads/"
            },
            { 
                name: "gsi software modules : samtools/1.9 bwa/0.7.12",
                url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
            },
            {   
                name: "gsi hg38 modules:  hg38-bwa-index/0.7.12",
                url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
            },
            {   
                name: "gsi hg19 modules:  hg19-bwa-index/0.7.12",
                url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
            },
            {   
                name: "gsi mm10 modules:  mm10-bwa-index/0.7.12",
                url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
            },
            {
                name: "bam-qc-metrics/0.2.5",
                url: "https://github.com/oicr-gsi/bam-qc-metrics.git"
            }

        ]
        output_meta: {
            deduplicatedBam: "Bam file after deduplication",
            statsEditDistance: "tsv file reports the (binned) average edit distance between the UMIs at each position",
            umiCountsPerPosition: "tsv file tabulates the counts for unique combinations of UMI and position",
            umiCounts: "tsv file provides UMI-level summary statistics",
            preDedupBamMetrics: "pre-collapse bamqc metrics",
            postDedupBamMetrics: "post-collapse bamqc metrics"
        }
    }

    output {
        File? deduplicatedBam = umiDedupedBam
        File statsEditDistance = statsMerge.statsEditDistance
        File umiCountsPerPosition = statsMerge.umiCountsPerPosition
        File umiCounts = statsMerge.umiCounts
        File? preDedupBamMetrics = preDedupBamQC.result
        File? postDedupBamMetrics = postDedupBamQC.result
    } 

    call getUMILengths {
      input:
          umiList = umiList    

    }

    Array[Int] umiLengths = getUMILengths.umiLengths

    scatter (fq in fastqInputs) {
        call extractUMIs { 
            input:
                umiList = umiList,
                outputPrefix = outputPrefix,
                fastq1 = fq.fastq1,
                fastq2 = fq.fastq2,
                pattern1 = pattern1,
                pattern2 = pattern2
            }
            call bwaMem.bwaMem {
                input:
                    fastqR1 = extractUMIs.fastqR1,
                    fastqR2 = extractUMIs.fastqR2,
                    readGroups = fq.readGroups,
                    outputFileNamePrefix = outputPrefix,
                    runBwaMem_bwaRef = resources[reference].bwaMem_runBwaMem_bwaRef,
                    runBwaMem_modules = resources[reference].bwaMem_runBwaMem_modules
            }
    }

    call bamMerge as mergeLibrary {
        input:
            outputPrefix = outputPrefix,
            Bams = bwaMem.bwaMemBam
    }

    if (doBamQC) {
        call bamQC.bamQC as preDedupBamQC {
            input:
                inputGroups = [{"bam": mergeLibrary.mergedBam, "bamIndex": mergeLibrary.mergedBai}],
                outputFileNamePrefix = "~{outputPrefix}.preDedup",
                mode = mode
          }
    }

    scatter (umiLength in umiLengths) {
        call bamSplitDeduplication {
            input:
                bamFile = mergeLibrary.mergedBam,
                umiLength = umiLength,
                outputPrefix = outputPrefix
        }
    }

    call bamMerge {
        input:
            outputPrefix = outputPrefix,
            Bams = bamSplitDeduplication.umiDedupBams
    }

    if (provisionBam) {
        File umiDedupedBam = bamMerge.mergedBam 
    }

    if (doBamQC) {
        call bamQC.bamQC as postDedupBamQC {
            input:
                inputGroups = [{"bam": mergeLibrary.mergedBam, "bamIndex": mergeLibrary.mergedBai}],
                outputFileNamePrefix = "~{outputPrefix}.postDedup",
                mode = mode
        }
    }

    call statsMerge {
        input:
            statsEditDistances = bamSplitDeduplication.dedupEditDistance,
            umiCountsPerPositions = bamSplitDeduplication.dedupUmiCountsPerPosition,
            umiCountsArray = bamSplitDeduplication.dedupUmiCounts,
            outputPrefix = outputPrefix
    }
}

task getUMILengths {
  input {
      File umiList
  }

  parameter_meta {
      umiList: "File with valid UMIs"
  }

  command <<<
      set -euo pipefail
      k=($(awk '{ match($1, "([ACTG])+"); print RLENGTH }' ~{umiList} | uniq))

      for i in ${k[@]}
      do
          for j in ${k[@]}
          do
              # adding 1 to account for period in UMI
              # 'ATG.ATCG'
              L+=($(($i+$j+1)))
          done
      done

      umiLengths=($(tr ' ' '\n' <<< "${L[@]}" | awk '!u[$0]++' | tr ' ' '\n'))
      printf "%s\n" "${umiLengths[@]}"

  >>>

  output {

    Array[Int] umiLengths = read_lines(stdout())

  }

  meta {
            output_meta: {
                umiLengths: "An integer array of UMI lengths"
            }
        }

}

task extractUMIs {
        input {
            String umiList
            String outputPrefix
            File fastq1
            File fastq2
            String pattern1
            String pattern2
            String modules = "barcodex-rs/0.1.2 rust/1.45.1"
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            umiList: "Reference file with valid UMIs"
            outputPrefix: "Specifies the start of the output files"
            fastq1: "FASTQ file containing read 1"
            fastq2: "FASTQ file containing read 2"
            pattern1: "UMI RegEx pattern 1"
            pattern2: "UMI RegEx pattern 2"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            set -euo pipefail
            barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
            --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
            --pattern2 '~{pattern2}' --r2-in ~{fastq2} 

            cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt

            tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
            echo "{$(sort -i tmp.txt)}" > new.txt
            tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            File fastqR1 = "~{outputPrefix}_R1.fastq.gz"
            File fastqR2 = "~{outputPrefix}_R2.fastq.gz"
            File discardR1 = "~{outputPrefix}_R1.discarded.fastq.gz"
            File discardR2 = "~{outputPrefix}_R2.discarded.fastq.gz"
            File extractR1 = "~{outputPrefix}_R1.extracted.fastq.gz"
            File extractR2 = "~{outputPrefix}_R2.extracted.fastq.gz"
            File umiCounts = "~{outputPrefix}_UMI_counts.json"
            File extractionMetrics = "~{outputPrefix}_extraction_metrics.json"
        }

        meta {
            output_meta: {
                fastqR1: "Read 1 fastq file with UMIs extracted",
                fastqR2: "Read 2 fastq file with UMIs extracted",
                discardR1: "Reads without a matching UMI pattern in read 1",
                discardR2: "Reads without a matching UMI pattern in read 2",
                extractR1: "Extracted reads (UMIs and any spacer sequences) from read 1",
                extractR2: "Extracted reads (UMIs and any spacer sequences) from read 2",
                umiCounts: "Record of UMI counts after extraction",
                extractionMetrics: "Metrics relating to extraction process"
            }
        }
}


task bamSplitDeduplication {
    input {
        Int umiLength
        File bamFile
        String modules = "umi-tools/1.0.0 samtools/1.9"
        String outputPrefix
        Int memory = 24
        Int timeout = 6
        String method = "directional"
        Int editDistanceThreshold = 1
    }

    parameter_meta {
        bamFile: "Bam file from bwaMem containing UMIs of varying lengths"
        umiLength: "Specifies the start of the output files"
        outputPrefix: "Specifies the start of the output files"
        modules: "Required environment modules"
        memory: "Memory allocated for this job"
        timeout: "Time in hours before task timeout"
        method: "What method to use to identify group of reads with the same (or similar) UMI(s)?"
        editDistanceThreshold: "Parametr for the adjacency and cluster methods, the threshold for the edit distance to connect two UMIs in the network."
    }

    command <<<
        set -euo pipefail
        samtools view -H ~{bamFile} > ~{outputPrefix}.~{umiLength}.sam
        samtools view ~{bamFile} | grep -P "^.*__\S{~{umiLength}}\t" >> ~{outputPrefix}.~{umiLength}.sam
        samtools view -Sb ~{outputPrefix}.~{umiLength}.sam > ~{outputPrefix}.~{umiLength}.bam

        samtools index ~{outputPrefix}.~{umiLength}.bam

        umi_tools dedup -I ~{outputPrefix}.~{umiLength}.bam \
        -S deduplicated.bam \
        --method=~{method} \
        --edit-distance-threshold=~{editDistanceThreshold} \
        --output-stats=deduplicated 
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File umiDedupBams = "deduplicated.bam"
        File dedupUmiCountsPerPosition = "deduplicated_per_umi_per_position.tsv"
        File dedupEditDistance = "deduplicated_edit_distance.tsv"
        File dedupUmiCounts = "deduplicated_per_umi.tsv"
    }

    meta {
        output_meta: {
            umiDedupBams: "Bam files with deduplicated UMIs of varying lengths",
            dedupUmiCountsPerPosition: "tsv file tabulates the counts for unique combinations of UMI and position"
        }
    }
}

task bamMerge {
    input {
        Array[File] Bams
        String modules = "samtools/1.9"
        Int memory = 24
        Int timeout = 6
        String outputPrefix
    }

    parameter_meta {
        Bams: "Input bam files"
        outputPrefix: "Prefix for output file"
        memory: "Memory allocated for indexing job"
        modules: "Required environment modules"
        timeout: "Hours before task timeout"
    }
    
    String resultMergedBam = "~{outputPrefix}.dedup.bam"

    command <<<        
        samtools merge -c ~{outputPrefix}.dedup.bam ~{sep=" " Bams}
        samtools index "~{outputPrefix}.dedup.bam"
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File mergedBam = "~{resultMergedBam}"
        File mergedBai = "~{outputPrefix}.dedup.bam.bai"
    }

    meta {
        output_meta: {
            mergedBam : "Deduplicated bam file"
        }
    }
}

task statsMerge {
    input{
        Array[File] umiCountsPerPositions
        Array[File] umiCountsArray
        Array[File] statsEditDistances
        Int memory = 16
        String outputPrefix
    }

    parameter_meta {
        umiCountsPerPositions: "An array of TSV files with umiCountsPerPosition"
        memory: "Memory allocated for this job"
    }

    command <<<
        set -euo pipefail
        statsEditDistances=(~{sep=" " statsEditDistances})
        umiCountsPerPositions=(~{sep=" " umiCountsPerPositions})
        umiCountsArray=(~{sep=" " umiCountsArray})

        length=${#statsEditDistances[@]}
        touch stats.tsv
        i=0
        while [ $i -lt $length ]
        do
            cat stats.tsv > tmp.tsv
            awk  'BEGIN {OFS='\t'} \
            {s1[$5]+=$1; s2[$5]+=$2; s3[$5]+=$3; s4[$5]+=$4} \
            END {for (k in s1) print s1[k]"\t"s2[k]"\t"s3[k]"\t"s4[k]"\t"k}' \
            tmp.tsv <(tail -n +2 ${statsEditDistances[i]}) > stats.tsv
            i=$(( $i+1 )) 
        done
        sort -k5 -o stats.tsv stats.tsv
        cat <(head -n 1 ${statsEditDistances[0]}) stats.tsv > ~{outputPrefix}.statsEditDistance.tsv

        length=${#umiCountsPerPositions[@]}
        touch umi.tsv
        i=0
        while [ $i -lt $length ]
        do
            cat umi.tsv > tmp.tsv
            awk  '{s2[$1]+=$2; s3[$1]+=$3} \
            END {for (k in s2) print k"\t"s2[k]"\t"s3[k]}' \
            tmp.tsv <(tail -n +2 ${umiCountsPerPositions[i]}) > umi.tsv
            i=$(( $i+1 ))
        done
        cat <(head -n 1 ${umiCountsPerPositions[0]}) umi.tsv > ~{outputPrefix}.umiCountsPerPosition.tsv

        length=${#umiCountsArray[@]}
        head -n 1 ${umiCountsArray[0]} > umiCounts.tsv
        i=0
        while [ $i -lt $length ]
        do
            cat umiCounts.tsv > tmp.tsv
            cat tmp.tsv <(tail -n +2 ${umiCountsArray[i]}) > ~{outputPrefix}.umiCounts.tsv
            i=$(( $i+1 ))
        done
    >>>

    runtime {
    memory: "~{memory}G"
    }

    output {
        File umiCountsPerPosition = "~{outputPrefix}.umiCountsPerPosition.tsv"
        File statsEditDistance = "~{outputPrefix}.statsEditDistance.tsv"
        File umiCounts= "~{outputPrefix}.umiCounts.tsv"
    }
}