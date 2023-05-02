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


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputPrefix`|String|"output"|Specifies the prefix of output files


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
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`mergeLibrary.modules`|String|"samtools/1.9"|Required environment modules
`mergeLibrary.memory`|Int|24|Memory allocated for indexing job
`mergeLibrary.timeout`|Int|6|Hours before task timeout
`bamSplitDeduplication.modules`|String|"umi-tools/1.0.0 samtools/1.9"|Required environment modules
`bamSplitDeduplication.memory`|Int|24|Memory allocated for this job
`bamSplitDeduplication.timeout`|Int|6|Time in hours before task timeout
`bamSplitDeduplication.method`|String|"directional"|What method to use to identify group of reads with the same (or similar) UMI(s)?
`bamSplitDeduplication.editDistanceThreshold`|Int|1|Parametr for the adjacency and cluster methods, the threshold for the edit distance to connect two UMIs in the network.
`bamMerge.modules`|String|"samtools/1.9"|Required environment modules
`bamMerge.memory`|Int|24|Memory allocated for indexing job
`bamMerge.timeout`|Int|6|Hours before task timeout
`statsMerge.memory`|Int|16|Memory allocated for this job


### Outputs

Output | Type | Description
---|---|---
`deduplicatedBam`|File|Bam file after deduplication
`statsEditDistance`|File|tsv file reports the (binned) average edit distance between the UMIs at each position
`umiCountsPerPosition`|File|tsv file tabulates the counts for unique combinations of UMI and position
`umiCounts`|File|tsv file provides UMI-level summary statistics


## Commands
   This section lists command(s) run by UmiCollapse workflow
   
   * Running UmiCollapse
   
  
   ```
   
               barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
               --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
               --pattern2 '~{pattern2}' --r2-in ~{fastq2} 
   
               cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt
   
               tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
               echo "{$(sort -i tmp.txt)}" > new.txt
               tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
  ```
   ```
           samtools view -H ~{bamFile} > ~{outputPrefix}.~{umiLength}.sam
           samtools view ~{bamFile} | grep -P "^.*__\S{~{umiLength}}\t" >> ~{outputPrefix}.~{umiLength}.sam
           samtools view $bamfile | grep -P "^.*__\S{~{7}\t" >> output.7.sam
           samtools view -Sb ~{outputPrefix}.~{umiLength}.sam > ~{outputPrefix}.~{umiLength}.bam
   
           samtools index ~{outputPrefix}.~{umiLength}.bam
   
           umi_tools dedup \
           -I ~{outputPrefix}.~{umiLength}.bam \
           -S deduplicated.bam \
           --output-stats=deduplicated \
           --log=deduplicated.log \
           --paired 
  ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
