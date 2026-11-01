# umiCollapse

## Overview

The incorporation of Unique Molecular Indices (UMIs) into sequenced reads allows for more accurate identification of PCR duplicates. This workflow extracts UMIs from the reads in fastq files, into the sequence identifier line. UMI identification is based on a known location and pattern and with the ability to match to a list of expected sequences. Reads are then aligned to a reference genome with bwa, and then the aligned sequence file is collapsed to remove duplicates with um-Tools

## Dependencies

* [barcodex-rs 0.1.2](https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz)
* [rust 1.2](https://www.rust-lang.org/tools/install)
* [umi-tools 1.1.1](https://github.com/CGATOxford/UMI-tools/archive/1.1.1.tar.gz)
* [samtools 1.9](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)
* [python 3.6](https://www.python.org/downloads/)
* [gsi software modules : samtools 1.9 bwa-mem2 2.2.1](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [bam-qc-metrics 0.2.8](https://github.com/oicr-gsi/bam-qc-metrics.git)
* [bwamem2 1.1.2](https://github.com/oicr-gsi/bwamem2)
* [bamQC 5.3.7](https://github.com/oicr-gsi/bamQC)


## Usage

### Cromwell
```
java -jar cromwell.jar run umiCollapse.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`umiList`|String|Reference file with valid UMIs
`fastqInputs`|Array[FastqInputs]|Array of fastq structs containing reads and readgroups
`pattern1`|String|UMI pattern 1
`pattern2`|String|UMI pattern 2
`reference`|String|Name and version of reference genome (hg19, hg38, hg38_noAlt; hg19_noAlt is not supported when doBamQC is true)
`preDedupBamQC.metadata`|Map[String,String]|JSON file containing metadata
`postDedupBamQC.metadata`|Map[String,String]|JSON file containing metadata


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputPrefix`|String|"output"|Specifies the prefix of output files
`doBamQC`|Boolean|false|Enable/disable bamQC process
`provisionBam`|Boolean|true|Enable/disable provision out umi deduplicated bam file


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractUMIs.modules`|String|"barcodex-rs/0.1.2 rust/1.45.1"|Required environment modules
`extractUMIs.memory`|Int|24|Memory allocated for this job
`extractUMIs.timeout`|Int|12|Time in hours before task timeout
`bwamem2.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwamem2.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwamem2.indexBam_timeout`|Int|48|Hours before task timeout
`bwamem2.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwamem2.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwamem2.bamMerge_timeout`|Int|72|Hours before task timeout
`bwamem2.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwamem2.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwamem2.runBwamem2_timeout`|Int|96|Hours before task timeout
`bwamem2.runBwamem2_jobMemory`|Int|32|Memory allocated for this job
`bwamem2.runBwamem2_threads`|Int|8|Requested CPU threads
`bwamem2.runBwamem2_addParam`|String?|None|Additional BWA parameters
`bwamem2.adapterTrimming_adapterTrim`|Boolean|true|If false, will not preform adapter trimming
`bwamem2.adapterTrimming_polyGTrim`|Int?|None|Number to pass to --nexseq-trim. Will address polyG trimming
`bwamem2.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwamem2.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwamem2.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwamem2.adapterTrimming_adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|Adapter sequence to trim from read 2
`bwamem2.adapterTrimming_adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|Adapter sequence to trim from read 1
`bwamem2.adapterTrimming_trimMinQuality`|Int|0|Minimum quality of read ends to keep
`bwamem2.adapterTrimming_trimMinLength`|Int|1|Minimum length of reads to keep
`bwamem2.adapterTrimming_umiLength`|Int|5|The number of bases to trim when doUMItrim is true. If the given length is positive, the bases are removed from the beginning of each read. If it is negative, the bases are removed from the end
`bwamem2.adapterTrimming_doUMItrim`|Boolean|false|If true, do umi trimming
`bwamem2.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwamem2.extractUMIs_timeout`|Int|12|Time in hours before task timeout
`bwamem2.extractUMIs_jobMemory`|Int|24|Memory allocated for this job
`bwamem2.extractUMIs_modules`|String|"barcodex-rs/0.1.2 rust/1.45.1"|Required environment modules
`bwamem2.extractUMIs_pattern2`|String|"(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"|UMI RegEx pattern 2
`bwamem2.extractUMIs_pattern1`|String|"(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"|UMI RegEx pattern 1
`bwamem2.extractUMIs_outputPrefix`|String|"extractUMIs_output"|Specifies the start of the output files
`bwamem2.extractUMIs_umiList`|String|"umiList"|Reference file with valid UMIs
`bwamem2.slicerR2_timeout`|Int|48|Hours before task timeout
`bwamem2.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwamem2.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwamem2.slicerR1_timeout`|Int|48|Hours before task timeout
`bwamem2.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwamem2.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwamem2.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwamem2.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwamem2.countChunkSize_modules`|String|"python/3.7"|Required environment modules
`bwamem2.numChunk`|Int|1|Number of chunks to split fastq file [1, no splitting]
`bwamem2.doUMIextract`|Boolean|false|If true, UMI will be extracted before alignment [false]
`bwamem2.doTrim`|Boolean|false|If true, adapters will be trimmed before alignment [false]
`bwamem2.numReads`|Int?|None|Number of reads
`mergeLibrary.modules`|String|"samtools/1.9"|Required environment modules
`mergeLibrary.memory`|Int|24|Memory allocated for indexing job
`mergeLibrary.timeout`|Int|6|Hours before task timeout
`preDedupBamQC.mergeReports_timeout`|Int|2|Timeout in hours for the merging task
`preDedupBamQC.mergeReports_jobMemory`|Int|4|RAM allocated to run the merging task
`preDedupBamQC.mergeReports_bamQCmerger`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_merger.py"|Path to the merger script
`preDedupBamQC.mergeReports_modules`|String|"bam-qc-metrics/0.2.8"|Runtime modules
`preDedupBamQC.mergedCoverageToHistogram_timeout`|Int|1|hours before task timeout
`preDedupBamQC.mergedCoverageToHistogram_threads`|Int|4|Requested CPU threads
`preDedupBamQC.mergedCoverageToHistogram_jobMemory`|Int|8|Memory allocated for this job
`preDedupBamQC.mergedCoverageToHistogram_coverageMerge`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_coverage_merger.py"|Path to coverage merging script
`preDedupBamQC.mergedCoverageToHistogram_modules`|String|"bam-qc-metrics/0.2.8"|required environment modules
`preDedupBamQC.markDuplicatesMerged_timeout`|Int|4|hours before task timeout
`preDedupBamQC.markDuplicatesMerged_threads`|Int|4|Requested CPU threads
`preDedupBamQC.markDuplicatesMerged_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.markDuplicatesMerged_modules`|String|"samblaster/0.1.26 samtools/1.16.1"|required environment modules
`preDedupBamQC.markDuplicatesMerged_additionalParameters`|String?|None|Any additional parameters for samblaster
`preDedupBamQC.bamQCMetrics_timeout`|Int|12|hours before task timeout
`preDedupBamQC.bamQCMetrics_jobMemory`|Int|8|Memory allocated for this job
`preDedupBamQC.bamQCMetrics_bamQClite`|String|"$BAM_QC_METRICS_ROOT/bin/run_bam_qc_lite.py"|Path to bamQC lite script
`preDedupBamQC.bamQCMetrics_modules`|String|"bam-qc-metrics/0.2.8"|required environment modules
`preDedupBamQC.bamQCMetrics_workflowVersion`|String|"5.3.5"|Workflow version to put into report
`preDedupBamQC.runWindowedMosdepth_timeout`|Int|4|hours before task timeout
`preDedupBamQC.runWindowedMosdepth_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.runWindowedMosdepth_targetBed`|String?|None|Optional target bed file
`preDedupBamQC.runWindowedMosdepth_additionalParameters`|String?|None|Additional parameters for coverage calculation
`preDedupBamQC.runWindowedMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`preDedupBamQC.downsampleBam_timeout`|Int|4|hours before task timeout
`preDedupBamQC.downsampleBam_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.downsampleBam_modules`|String|"samtools/1.16.1"|required environment modules
`preDedupBamQC.getUniqueReadsCount_timeout`|Int|18|hours before task timeout
`preDedupBamQC.getUniqueReadsCount_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.getUniqueReadsCount_modules`|String|"samtools/1.16.1"|required environment modules
`preDedupBamQC.markDuplicates_timeout`|Int|4|hours before task timeout
`preDedupBamQC.markDuplicates_threads`|Int|4|Requested CPU threads
`preDedupBamQC.markDuplicates_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.markDuplicates_modules`|String|"samblaster/0.1.26 samtools/1.16.1"|required environment modules
`preDedupBamQC.markDuplicates_additionalParameters`|String?|None|Any additional parameters for samblaster
`preDedupBamQC.cumulativeDistToHistogram_timeout`|Int|1|hours before task timeout
`preDedupBamQC.cumulativeDistToHistogram_threads`|Int|4|Requested CPU threads
`preDedupBamQC.cumulativeDistToHistogram_jobMemory`|Int|8|Memory allocated for this job
`preDedupBamQC.cumulativeDistToHistogram_outFileName`|String|"coverage_histogram.json"|Output file name, default coverage_histogram.json
`preDedupBamQC.cumulativeDistToHistogram_coverageHistogram`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_coverage_histogram.py"|Path to script generating coverage histogram
`preDedupBamQC.cumulativeDistToHistogram_modules`|String|"bam-qc-metrics/0.2.8"|required environment modules
`preDedupBamQC.runBedtoolsIntersect_timeout`|Int|12|hours before task timeout
`preDedupBamQC.runBedtoolsIntersect_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.runBedtoolsIntersect_modules`|String|"samtools/1.16.1 bedtools/2.27"|required environment modules
`preDedupBamQC.runMosdepth_timeout`|Int|4|hours before task timeout
`preDedupBamQC.runMosdepth_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.runMosdepth_additionalParameters`|String?|None|Additional parameters for coverage calculation
`preDedupBamQC.runMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`preDedupBamQC.samstats_timeout`|Int|24|hours before task timeout
`preDedupBamQC.samstats_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.samstats_additionalParameters`|String?|None|Additional parameters for samtools, usually filtering flag
`preDedupBamQC.targetBed`|String?|None|Path to optional target bed file
`preDedupBamQC.filterTargetedCoverage`|Boolean|true|Flag for enabling filtering when calculating reads on target
`preDedupBamQC.downsampleToReads`|Int|500000|Downsample to this many reads when running unique read count, duplicate rate calculation and CIGAR analysis
`preDedupBamQC.coverageWindow`|Int|1000|Coverage window to use with mosdepth for making coverage histogram, default is 1000b
`bamSplitDeduplication.modules`|String|"umi-tools/1.0.0 samtools/1.9"|Required environment modules
`bamSplitDeduplication.memory`|Int|24|Memory allocated for this job
`bamSplitDeduplication.timeout`|Int|6|Time in hours before task timeout
`bamSplitDeduplication.method`|String|"directional"|What method to use to identify group of reads with the same (or similar) UMI(s)?
`bamSplitDeduplication.editDistanceThreshold`|Int|1|Parametr for the adjacency and cluster methods, the threshold for the edit distance to connect two UMIs in the network.
`bamMerge.modules`|String|"samtools/1.9"|Required environment modules
`bamMerge.memory`|Int|24|Memory allocated for indexing job
`bamMerge.timeout`|Int|6|Hours before task timeout
`postDedupBamQC.mergeReports_timeout`|Int|2|Timeout in hours for the merging task
`postDedupBamQC.mergeReports_jobMemory`|Int|4|RAM allocated to run the merging task
`postDedupBamQC.mergeReports_bamQCmerger`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_merger.py"|Path to the merger script
`postDedupBamQC.mergeReports_modules`|String|"bam-qc-metrics/0.2.8"|Runtime modules
`postDedupBamQC.mergedCoverageToHistogram_timeout`|Int|1|hours before task timeout
`postDedupBamQC.mergedCoverageToHistogram_threads`|Int|4|Requested CPU threads
`postDedupBamQC.mergedCoverageToHistogram_jobMemory`|Int|8|Memory allocated for this job
`postDedupBamQC.mergedCoverageToHistogram_coverageMerge`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_coverage_merger.py"|Path to coverage merging script
`postDedupBamQC.mergedCoverageToHistogram_modules`|String|"bam-qc-metrics/0.2.8"|required environment modules
`postDedupBamQC.markDuplicatesMerged_timeout`|Int|4|hours before task timeout
`postDedupBamQC.markDuplicatesMerged_threads`|Int|4|Requested CPU threads
`postDedupBamQC.markDuplicatesMerged_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.markDuplicatesMerged_modules`|String|"samblaster/0.1.26 samtools/1.16.1"|required environment modules
`postDedupBamQC.markDuplicatesMerged_additionalParameters`|String?|None|Any additional parameters for samblaster
`postDedupBamQC.bamQCMetrics_timeout`|Int|12|hours before task timeout
`postDedupBamQC.bamQCMetrics_jobMemory`|Int|8|Memory allocated for this job
`postDedupBamQC.bamQCMetrics_bamQClite`|String|"$BAM_QC_METRICS_ROOT/bin/run_bam_qc_lite.py"|Path to bamQC lite script
`postDedupBamQC.bamQCMetrics_modules`|String|"bam-qc-metrics/0.2.8"|required environment modules
`postDedupBamQC.bamQCMetrics_workflowVersion`|String|"5.3.5"|Workflow version to put into report
`postDedupBamQC.runWindowedMosdepth_timeout`|Int|4|hours before task timeout
`postDedupBamQC.runWindowedMosdepth_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.runWindowedMosdepth_targetBed`|String?|None|Optional target bed file
`postDedupBamQC.runWindowedMosdepth_additionalParameters`|String?|None|Additional parameters for coverage calculation
`postDedupBamQC.runWindowedMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`postDedupBamQC.downsampleBam_timeout`|Int|4|hours before task timeout
`postDedupBamQC.downsampleBam_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.downsampleBam_modules`|String|"samtools/1.16.1"|required environment modules
`postDedupBamQC.getUniqueReadsCount_timeout`|Int|18|hours before task timeout
`postDedupBamQC.getUniqueReadsCount_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.getUniqueReadsCount_modules`|String|"samtools/1.16.1"|required environment modules
`postDedupBamQC.markDuplicates_timeout`|Int|4|hours before task timeout
`postDedupBamQC.markDuplicates_threads`|Int|4|Requested CPU threads
`postDedupBamQC.markDuplicates_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.markDuplicates_modules`|String|"samblaster/0.1.26 samtools/1.16.1"|required environment modules
`postDedupBamQC.markDuplicates_additionalParameters`|String?|None|Any additional parameters for samblaster
`postDedupBamQC.cumulativeDistToHistogram_timeout`|Int|1|hours before task timeout
`postDedupBamQC.cumulativeDistToHistogram_threads`|Int|4|Requested CPU threads
`postDedupBamQC.cumulativeDistToHistogram_jobMemory`|Int|8|Memory allocated for this job
`postDedupBamQC.cumulativeDistToHistogram_outFileName`|String|"coverage_histogram.json"|Output file name, default coverage_histogram.json
`postDedupBamQC.cumulativeDistToHistogram_coverageHistogram`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_coverage_histogram.py"|Path to script generating coverage histogram
`postDedupBamQC.cumulativeDistToHistogram_modules`|String|"bam-qc-metrics/0.2.8"|required environment modules
`postDedupBamQC.runBedtoolsIntersect_timeout`|Int|12|hours before task timeout
`postDedupBamQC.runBedtoolsIntersect_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.runBedtoolsIntersect_modules`|String|"samtools/1.16.1 bedtools/2.27"|required environment modules
`postDedupBamQC.runMosdepth_timeout`|Int|4|hours before task timeout
`postDedupBamQC.runMosdepth_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.runMosdepth_additionalParameters`|String?|None|Additional parameters for coverage calculation
`postDedupBamQC.runMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`postDedupBamQC.samstats_timeout`|Int|24|hours before task timeout
`postDedupBamQC.samstats_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.samstats_additionalParameters`|String?|None|Additional parameters for samtools, usually filtering flag
`postDedupBamQC.targetBed`|String?|None|Path to optional target bed file
`postDedupBamQC.filterTargetedCoverage`|Boolean|true|Flag for enabling filtering when calculating reads on target
`postDedupBamQC.downsampleToReads`|Int|500000|Downsample to this many reads when running unique read count, duplicate rate calculation and CIGAR analysis
`postDedupBamQC.coverageWindow`|Int|1000|Coverage window to use with mosdepth for making coverage histogram, default is 1000b
`statsMerge.memory`|Int|16|Memory allocated for this job


### Outputs

Output | Type | Description | Labels
---|---|---|---
`deduplicatedBam`|File?|Bam file after deduplication|vidarr_label: deduplicatedBam
`statsEditDistance`|File|tsv file reports the (binned) average edit distance between the UMIs at each position|vidarr_label: statsEditDistance
`umiCountsPerPosition`|File|tsv file tabulates the counts for unique combinations of UMI and position|vidarr_label: umiCountsPerPosition
`umiCounts`|File|tsv file provides UMI-level summary statistics|vidarr_label: umiCounts
`preDedupBamMetrics`|File?|pre-collapse bamqc metrics|vidarr_label: preDedupBamMetrics
`postDedupBamMetrics`|File?|post-collapse bamqc metrics|vidarr_label: postDedupBamMetrics


## Commands
This section lists command(s) run by umiCollapse workflow

* Running umiCollapse

```
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

      umiLengths=($(tr ' ' '\n' ``` "${L[@]}" | awk '!u[$0]++' | tr ' ' '\n'))
      printf "%s\n" "${umiLengths[@]}"

```
```
            set -euo pipefail
            barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
            --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
            --pattern2 '~{pattern2}' --r2-in ~{fastq2} 

            cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt

            tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
            echo "{$(sort -i tmp.txt)}" > new.txt
            tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
```
```
        set -euo pipefail
        samtools view -H ~{bamFile} > ~{outputPrefix}.~{umiLength}.sam
        samtools view -F 2304 ~{bamFile} | grep -P "^.*__\S{~{umiLength}}\t" >> ~{outputPrefix}.~{umiLength}.sam
        samtools view -Sb ~{outputPrefix}.~{umiLength}.sam > ~{outputPrefix}.~{umiLength}.bam

        samtools index ~{outputPrefix}.~{umiLength}.bam

        umi_tools dedup --paired -I ~{outputPrefix}.~{umiLength}.bam \
        -S deduplicated.bam \
        --method=~{method} \
        --edit-distance-threshold=~{editDistanceThreshold} \
        --output-stats=deduplicated 
```
```        
        samtools merge -c ~{outputPrefix}.dedup.bam ~{sep=" " Bams}
        samtools index "~{outputPrefix}.dedup.bam"
```
```
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
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
