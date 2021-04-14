# Intro

The Nextflow script runs SGA a list of sample reads.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimited text file and should be specified in `nextflow.config`.  For the SGA run, SampleID, FastqR1 and FastqR2 columns are required.

- FastqR1 and FastqR2 should contain the full path to the sample Fastq reads.

| SampleID | FastqR1 | FastqR2 |
| -------- | ------- | ------- |
| A01      | /path-to/A01_R1.fastq.gz       | /path-to/A01_R2.fastq.gz       |

## Examples

Attached is `NA12878.samplesheet.tsv` an example sheet and `nextflow.config` to show configuration settings.

## To run

For each dataset
1) Create your sample sheet. E.g. `NA12878.samplesheet.tsv`.
2) Modify your `nextflow.config` to read the `NA12878.samplesheet.tsv` and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/1kg/datasets/NA12878/nextflow-out"`
3) Run the workflow
```
nextflow -log nextflow.log run -w /cbio/projects/012/stream1/team/gerrit/assembly/work -c /cbio/projects/012/stream1/team/gerrit/assembly/dec-2020-hackathon-stream1/assembly/nextflow.config /cbio/projects/012/stream1/team/gerrit/assembly/dec-2020-hackathon-stream1/assembly/main.nf -with-report report.html -with-timeline timeline.html -profile ilifu -resume
```

## Output

The output directory will contain per sample directories. Each sample directory will contain

1. `NA12878.*` - All SGA created files
2. `NA12878.assemble`, `NA12878.assemble-graph.asqg.gz`, `NA12878.assemble-variants.fa`, `NA12878.correct.filter.pass.merged.rmdup.asqg.gz`  - output files generated from HUPAN pipeline
