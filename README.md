# umiCollapse

The incorporation of Unique Molecular Indices (UMIs) into sequenced reads allows for more accurate identification of PCR duplicates. This workflow extracts UMIs from the reads in fastq files, into the sequence identifier line. UMI identification is based on a known location and pattern and with the ability to match to a list of expected sequences. Reads are then aligned to a reference genome with bwa, and then the aligned sequence file is collapsed to remove duplicates with um-Tools

## Overview

## Dependencies

* [barcodex-rs 0.1.2](https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz)
* [rust 1.2](https://www.rust-lang.org/tools/install)
* [umi-tools 1.1.1](https://github.com/CGATOxford/UMI-tools/archive/1.1.1.tar.gz)
* [bwa 0.7.12](https://github.com/lh3/bwa/archive/0.7.12.tar.gz)
* [samtools 1.9](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)
* [python 3.6](https://www.python.org/downloads/)
* [gsi software modules : samtools 1.9 bwa 0.7.12](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi hg38 modules:  hg38-bwa-index 0.7.12](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi hg19 modules:  hg19-bwa-index 0.7.12](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi mm10 modules:  mm10-bwa-index 0.7.12](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [bam-qc-metrics 0.2.5](https://github.com/oicr-gsi/bam-qc-metrics.git)
* [bwaMem 2.2.1](https://github.com/oicr-gsi/bwaMem)
* [bamQC 5.1.3](https://github.com/oicr-gsi/bamQC)

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
`reference`|String|Name and version of reference genome
`mode`|String|running mode for the workflow, only allow value 'lane_level' and 'call_ready'
`bwaMem.reference`|String|The genome reference build. For example: hg19, hg38, mm10
`preDedupBamQC.bamQCMetrics_workflowVersion`|String|Workflow version string
`preDedupBamQC.bamQCMetrics_refSizesBed`|String|Path to human genome BED reference with chromosome sizes
`preDedupBamQC.bamQCMetrics_refFasta`|String|Path to human genome FASTA reference
`preDedupBamQC.metadata`|Map[String,String]|JSON file containing metadata
`postDedupBamQC.bamQCMetrics_workflowVersion`|String|Workflow version string
`postDedupBamQC.bamQCMetrics_refSizesBed`|String|Path to human genome BED reference with chromosome sizes
`postDedupBamQC.bamQCMetrics_refFasta`|String|Path to human genome FASTA reference
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
`bwaMem.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.indexBam_timeout`|Int|48|Hours before task timeout
`bwaMem.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwaMem.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.bamMerge_timeout`|Int|72|Hours before task timeout
`bwaMem.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwaMem.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwaMem.runBwaMem_timeout`|Int|96|Hours before task timeout
`bwaMem.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`bwaMem.runBwaMem_threads`|Int|8|Requested CPU threads
`bwaMem.runBwaMem_addParam`|String?|None|Additional BWA parameters
`bwaMem.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwaMem.adapterTrimming_adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|Adapter sequence to trim from read 2
`bwaMem.adapterTrimming_adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|Adapter sequence to trim from read 1
`bwaMem.adapterTrimming_trimMinQuality`|Int|0|Minimum quality of read ends to keep
`bwaMem.adapterTrimming_trimMinLength`|Int|1|Minimum length of reads to keep
`bwaMem.adapterTrimming_umiLength`|Int|5|The number of bases to trim when doUMItrim is true. If the given length is positive, the bases are removed from the beginning of each read. If it is negative, the bases are removed from the end
`bwaMem.adapterTrimming_doUMItrim`|Boolean|false|If true, do umi trimming
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.extractUMIs_timeout`|Int|12|Time in hours before task timeout
`bwaMem.extractUMIs_jobMemory`|Int|24|Memory allocated for this job
`bwaMem.extractUMIs_modules`|String|"barcodex-rs/0.1.2 rust/1.45.1"|Required environment modules
`bwaMem.extractUMIs_pattern2`|String|"(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"|UMI RegEx pattern 2
`bwaMem.extractUMIs_pattern1`|String|"(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"|UMI RegEx pattern 1
`bwaMem.extractUMIs_outputPrefix`|String|"extractUMIs_output"|Specifies the start of the output files
`bwaMem.extractUMIs_umiList`|String|"umiList"|Reference file with valid UMIs
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.countChunkSize_modules`|String|"python/3.7"|Required environment modules
`bwaMem.numChunk`|Int|1|Number of chunks to split fastq file [1, no splitting]
`bwaMem.doUMIextract`|Boolean|false|If true, UMI will be extracted before alignment [false]
`bwaMem.doTrim`|Boolean|false|If true, adapters will be trimmed before alignment [false]
`bwaMem.numReads`|Int?|None|Number of reads
`mergeLibrary.modules`|String|"samtools/1.9"|Required environment modules
`mergeLibrary.memory`|Int|24|Memory allocated for indexing job
`mergeLibrary.timeout`|Int|6|Hours before task timeout
`preDedupBamQC.collateResults_timeout`|Int|1|hours before task timeout
`preDedupBamQC.collateResults_threads`|Int|4|Requested CPU threads
`preDedupBamQC.collateResults_jobMemory`|Int|8|Memory allocated for this job
`preDedupBamQC.collateResults_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.cumulativeDistToHistogram_timeout`|Int|1|hours before task timeout
`preDedupBamQC.cumulativeDistToHistogram_threads`|Int|4|Requested CPU threads
`preDedupBamQC.cumulativeDistToHistogram_jobMemory`|Int|8|Memory allocated for this job
`preDedupBamQC.cumulativeDistToHistogram_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.runMosdepth_timeout`|Int|4|hours before task timeout
`preDedupBamQC.runMosdepth_threads`|Int|4|Requested CPU threads
`preDedupBamQC.runMosdepth_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.runMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`preDedupBamQC.bamQCMetrics_timeout`|Int|4|hours before task timeout
`preDedupBamQC.bamQCMetrics_threads`|Int|4|Requested CPU threads
`preDedupBamQC.bamQCMetrics_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.bamQCMetrics_modules`|String|"bam-qc-metrics/0.2.5"|required environment modules
`preDedupBamQC.bamQCMetrics_normalInsertMax`|Int|1500|Maximum of expected insert size range
`preDedupBamQC.markDuplicates_timeout`|Int|4|hours before task timeout
`preDedupBamQC.markDuplicates_threads`|Int|4|Requested CPU threads
`preDedupBamQC.markDuplicates_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.markDuplicates_modules`|String|"picard/2.21.2"|required environment modules
`preDedupBamQC.markDuplicates_picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR
`preDedupBamQC.markDuplicates_opticalDuplicatePixelDistance`|Int|100|Maximum offset between optical duplicate clusters
`preDedupBamQC.downsampleRegion_timeout`|Int|4|hours before task timeout
`preDedupBamQC.downsampleRegion_threads`|Int|4|Requested CPU threads
`preDedupBamQC.downsampleRegion_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.downsampleRegion_modules`|String|"samtools/1.14"|required environment modules
`preDedupBamQC.downsample_timeout`|Int|4|hours before task timeout
`preDedupBamQC.downsample_threads`|Int|4|Requested CPU threads
`preDedupBamQC.downsample_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.downsample_modules`|String|"samtools/1.14"|required environment modules
`preDedupBamQC.downsample_randomSeed`|Int|42|Random seed for pre-downsampling (if any)
`preDedupBamQC.downsample_downsampleSuffix`|String|"downsampled.bam"|Suffix for output file
`preDedupBamQC.findDownsampleParamsMarkDup_timeout`|Int|4|hours before task timeout
`preDedupBamQC.findDownsampleParamsMarkDup_threads`|Int|4|Requested CPU threads
`preDedupBamQC.findDownsampleParamsMarkDup_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.findDownsampleParamsMarkDup_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.findDownsampleParamsMarkDup_customRegions`|String|""|Custom downsample regions; overrides chromosome and interval parameters
`preDedupBamQC.findDownsampleParamsMarkDup_intervalStart`|Int|100000|Start of interval in each chromosome, for very large BAMs
`preDedupBamQC.findDownsampleParamsMarkDup_baseInterval`|Int|15000|Base width of interval in each chromosome, for very large BAMs
`preDedupBamQC.findDownsampleParamsMarkDup_chromosomes`|Array[String]|["chr12", "chr13", "chrXII", "chrXIII"]|Array of chromosome identifiers for downsampled subset
`preDedupBamQC.findDownsampleParamsMarkDup_threshold`|Int|10000000|Minimum number of reads to conduct downsampling
`preDedupBamQC.findDownsampleParams_timeout`|Int|4|hours before task timeout
`preDedupBamQC.findDownsampleParams_threads`|Int|4|Requested CPU threads
`preDedupBamQC.findDownsampleParams_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.findDownsampleParams_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.findDownsampleParams_preDSMultiplier`|Float|1.5|Determines target size for pre-downsampled set (if any). Must have (preDSMultiplier) < (minReadsRelative).
`preDedupBamQC.findDownsampleParams_precision`|Int|8|Number of decimal places in fraction for pre-downsampling
`preDedupBamQC.findDownsampleParams_minReadsRelative`|Int|2|Minimum value of (inputReads)/(targetReads) to allow pre-downsampling
`preDedupBamQC.findDownsampleParams_minReadsAbsolute`|Int|10000|Minimum value of targetReads to allow pre-downsampling
`preDedupBamQC.findDownsampleParams_targetReads`|Int|100000|Desired number of reads in downsampled output
`preDedupBamQC.updateMetadata_timeout`|Int|4|hours before task timeout
`preDedupBamQC.updateMetadata_threads`|Int|4|Requested CPU threads
`preDedupBamQC.updateMetadata_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.updateMetadata_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.filter_timeout`|Int|4|hours before task timeout
`preDedupBamQC.filter_threads`|Int|4|Requested CPU threads
`preDedupBamQC.filter_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.filter_modules`|String|"samtools/1.14"|required environment modules
`preDedupBamQC.filter_minQuality`|Int|30|Minimum alignment quality to pass filter
`preDedupBamQC.mergeFiles_modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`preDedupBamQC.mergeFiles_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`preDedupBamQC.mergeFiles_cores`|Int|1|The number of cores to allocate to the job.
`preDedupBamQC.mergeFiles_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`preDedupBamQC.mergeFiles_jobMemory`|Int|24|Memory allocated to job (in GB).
`preDedupBamQC.mergeFiles_suffix`|String|".merge"|suffix to use for merged bam
`preDedupBamQC.mergeSplitByIntervalFiles_modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`preDedupBamQC.mergeSplitByIntervalFiles_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`preDedupBamQC.mergeSplitByIntervalFiles_cores`|Int|1|The number of cores to allocate to the job.
`preDedupBamQC.mergeSplitByIntervalFiles_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`preDedupBamQC.mergeSplitByIntervalFiles_jobMemory`|Int|24|Memory allocated to job (in GB).
`preDedupBamQC.mergeSplitByIntervalFiles_suffix`|String|".merge"|suffix to use for merged bam
`preDedupBamQC.preFilter_timeout`|Int|4|hours before task timeout
`preDedupBamQC.preFilter_threads`|Int|4|Requested CPU threads
`preDedupBamQC.preFilter_minMemory`|Int|2|Minimum amount of RAM allocated to the task
`preDedupBamQC.preFilter_jobMemory`|Int|6|Memory allocated for this job
`preDedupBamQC.preFilter_modules`|String|"samtools/1.14"|required environment modules
`preDedupBamQC.preFilter_filterAdditionalParams`|String?|None|Additional parameters to pass to samtools.
`preDedupBamQC.preFilter_minMapQuality`|Int?|None|Minimum alignment quality to pass filter
`preDedupBamQC.preFilter_filterFlags`|Int|260|Samtools filter flags to apply.
`preDedupBamQC.getChrCoefficient_modules`|String|"samtools/1.14"|Names and versions of modules to load
`preDedupBamQC.getChrCoefficient_timeout`|Int|1|Hours before task timeout
`preDedupBamQC.getChrCoefficient_memory`|Int|2|Memory allocated for this job
`preDedupBamQC.splitStringToArray_modules`|String|""|Environment module name and version to load (space separated) before command execution.
`preDedupBamQC.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`preDedupBamQC.splitStringToArray_cores`|Int|1|The number of cores to allocate to the job.
`preDedupBamQC.splitStringToArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`preDedupBamQC.splitStringToArray_recordSeparator`|String|"+"|Record separator - a delimiter for joining records
`preDedupBamQC.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`preDedupBamQC.intervalsToParallelizeByString`|String|"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"|Comma separated list of intervals to split by (e.g. chr1,chr2,chr3,chr4).
`bamSplitDeduplication.modules`|String|"umi-tools/1.0.0 samtools/1.9"|Required environment modules
`bamSplitDeduplication.memory`|Int|24|Memory allocated for this job
`bamSplitDeduplication.timeout`|Int|6|Time in hours before task timeout
`bamSplitDeduplication.method`|String|"directional"|What method to use to identify group of reads with the same (or similar) UMI(s)?
`bamSplitDeduplication.editDistanceThreshold`|Int|1|Parametr for the adjacency and cluster methods, the threshold for the edit distance to connect two UMIs in the network.
`bamMerge.modules`|String|"samtools/1.9"|Required environment modules
`bamMerge.memory`|Int|24|Memory allocated for indexing job
`bamMerge.timeout`|Int|6|Hours before task timeout
`postDedupBamQC.collateResults_timeout`|Int|1|hours before task timeout
`postDedupBamQC.collateResults_threads`|Int|4|Requested CPU threads
`postDedupBamQC.collateResults_jobMemory`|Int|8|Memory allocated for this job
`postDedupBamQC.collateResults_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.cumulativeDistToHistogram_timeout`|Int|1|hours before task timeout
`postDedupBamQC.cumulativeDistToHistogram_threads`|Int|4|Requested CPU threads
`postDedupBamQC.cumulativeDistToHistogram_jobMemory`|Int|8|Memory allocated for this job
`postDedupBamQC.cumulativeDistToHistogram_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.runMosdepth_timeout`|Int|4|hours before task timeout
`postDedupBamQC.runMosdepth_threads`|Int|4|Requested CPU threads
`postDedupBamQC.runMosdepth_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.runMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`postDedupBamQC.bamQCMetrics_timeout`|Int|4|hours before task timeout
`postDedupBamQC.bamQCMetrics_threads`|Int|4|Requested CPU threads
`postDedupBamQC.bamQCMetrics_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.bamQCMetrics_modules`|String|"bam-qc-metrics/0.2.5"|required environment modules
`postDedupBamQC.bamQCMetrics_normalInsertMax`|Int|1500|Maximum of expected insert size range
`postDedupBamQC.markDuplicates_timeout`|Int|4|hours before task timeout
`postDedupBamQC.markDuplicates_threads`|Int|4|Requested CPU threads
`postDedupBamQC.markDuplicates_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.markDuplicates_modules`|String|"picard/2.21.2"|required environment modules
`postDedupBamQC.markDuplicates_picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR
`postDedupBamQC.markDuplicates_opticalDuplicatePixelDistance`|Int|100|Maximum offset between optical duplicate clusters
`postDedupBamQC.downsampleRegion_timeout`|Int|4|hours before task timeout
`postDedupBamQC.downsampleRegion_threads`|Int|4|Requested CPU threads
`postDedupBamQC.downsampleRegion_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.downsampleRegion_modules`|String|"samtools/1.14"|required environment modules
`postDedupBamQC.downsample_timeout`|Int|4|hours before task timeout
`postDedupBamQC.downsample_threads`|Int|4|Requested CPU threads
`postDedupBamQC.downsample_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.downsample_modules`|String|"samtools/1.14"|required environment modules
`postDedupBamQC.downsample_randomSeed`|Int|42|Random seed for pre-downsampling (if any)
`postDedupBamQC.downsample_downsampleSuffix`|String|"downsampled.bam"|Suffix for output file
`postDedupBamQC.findDownsampleParamsMarkDup_timeout`|Int|4|hours before task timeout
`postDedupBamQC.findDownsampleParamsMarkDup_threads`|Int|4|Requested CPU threads
`postDedupBamQC.findDownsampleParamsMarkDup_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.findDownsampleParamsMarkDup_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.findDownsampleParamsMarkDup_customRegions`|String|""|Custom downsample regions; overrides chromosome and interval parameters
`postDedupBamQC.findDownsampleParamsMarkDup_intervalStart`|Int|100000|Start of interval in each chromosome, for very large BAMs
`postDedupBamQC.findDownsampleParamsMarkDup_baseInterval`|Int|15000|Base width of interval in each chromosome, for very large BAMs
`postDedupBamQC.findDownsampleParamsMarkDup_chromosomes`|Array[String]|["chr12", "chr13", "chrXII", "chrXIII"]|Array of chromosome identifiers for downsampled subset
`postDedupBamQC.findDownsampleParamsMarkDup_threshold`|Int|10000000|Minimum number of reads to conduct downsampling
`postDedupBamQC.findDownsampleParams_timeout`|Int|4|hours before task timeout
`postDedupBamQC.findDownsampleParams_threads`|Int|4|Requested CPU threads
`postDedupBamQC.findDownsampleParams_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.findDownsampleParams_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.findDownsampleParams_preDSMultiplier`|Float|1.5|Determines target size for pre-downsampled set (if any). Must have (preDSMultiplier) < (minReadsRelative).
`postDedupBamQC.findDownsampleParams_precision`|Int|8|Number of decimal places in fraction for pre-downsampling
`postDedupBamQC.findDownsampleParams_minReadsRelative`|Int|2|Minimum value of (inputReads)/(targetReads) to allow pre-downsampling
`postDedupBamQC.findDownsampleParams_minReadsAbsolute`|Int|10000|Minimum value of targetReads to allow pre-downsampling
`postDedupBamQC.findDownsampleParams_targetReads`|Int|100000|Desired number of reads in downsampled output
`postDedupBamQC.updateMetadata_timeout`|Int|4|hours before task timeout
`postDedupBamQC.updateMetadata_threads`|Int|4|Requested CPU threads
`postDedupBamQC.updateMetadata_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.updateMetadata_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.filter_timeout`|Int|4|hours before task timeout
`postDedupBamQC.filter_threads`|Int|4|Requested CPU threads
`postDedupBamQC.filter_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.filter_modules`|String|"samtools/1.14"|required environment modules
`postDedupBamQC.filter_minQuality`|Int|30|Minimum alignment quality to pass filter
`postDedupBamQC.mergeFiles_modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`postDedupBamQC.mergeFiles_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`postDedupBamQC.mergeFiles_cores`|Int|1|The number of cores to allocate to the job.
`postDedupBamQC.mergeFiles_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`postDedupBamQC.mergeFiles_jobMemory`|Int|24|Memory allocated to job (in GB).
`postDedupBamQC.mergeFiles_suffix`|String|".merge"|suffix to use for merged bam
`postDedupBamQC.mergeSplitByIntervalFiles_modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`postDedupBamQC.mergeSplitByIntervalFiles_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`postDedupBamQC.mergeSplitByIntervalFiles_cores`|Int|1|The number of cores to allocate to the job.
`postDedupBamQC.mergeSplitByIntervalFiles_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`postDedupBamQC.mergeSplitByIntervalFiles_jobMemory`|Int|24|Memory allocated to job (in GB).
`postDedupBamQC.mergeSplitByIntervalFiles_suffix`|String|".merge"|suffix to use for merged bam
`postDedupBamQC.preFilter_timeout`|Int|4|hours before task timeout
`postDedupBamQC.preFilter_threads`|Int|4|Requested CPU threads
`postDedupBamQC.preFilter_minMemory`|Int|2|Minimum amount of RAM allocated to the task
`postDedupBamQC.preFilter_jobMemory`|Int|6|Memory allocated for this job
`postDedupBamQC.preFilter_modules`|String|"samtools/1.14"|required environment modules
`postDedupBamQC.preFilter_filterAdditionalParams`|String?|None|Additional parameters to pass to samtools.
`postDedupBamQC.preFilter_minMapQuality`|Int?|None|Minimum alignment quality to pass filter
`postDedupBamQC.preFilter_filterFlags`|Int|260|Samtools filter flags to apply.
`postDedupBamQC.getChrCoefficient_modules`|String|"samtools/1.14"|Names and versions of modules to load
`postDedupBamQC.getChrCoefficient_timeout`|Int|1|Hours before task timeout
`postDedupBamQC.getChrCoefficient_memory`|Int|2|Memory allocated for this job
`postDedupBamQC.splitStringToArray_modules`|String|""|Environment module name and version to load (space separated) before command execution.
`postDedupBamQC.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`postDedupBamQC.splitStringToArray_cores`|Int|1|The number of cores to allocate to the job.
`postDedupBamQC.splitStringToArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`postDedupBamQC.splitStringToArray_recordSeparator`|String|"+"|Record separator - a delimiter for joining records
`postDedupBamQC.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`postDedupBamQC.intervalsToParallelizeByString`|String|"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"|Comma separated list of intervals to split by (e.g. chr1,chr2,chr3,chr4).
`statsMerge.memory`|Int|16|Memory allocated for this job


### Outputs

Output | Type | Description
---|---|---
`deduplicatedBam`|File?|Bam file after deduplication
`statsEditDistance`|File|tsv file reports the (binned) average edit distance between the UMIs at each position
`umiCountsPerPosition`|File|tsv file tabulates the counts for unique combinations of UMI and position
`umiCounts`|File|tsv file provides UMI-level summary statistics
`preDedupBamMetrics`|File?|pre-collapse bamqc metrics
`postDedupBamMetrics`|File?|post-collapse bamqc metrics

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
