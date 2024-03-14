#### *Gamaliel López Leal*

#### Centro de Investigación en Dinámica Celular

***Mail:*** gamaliel.lopez@uaem.edu.mx

***Research Gate:***
*https://www.researchgate.net/profile/Gamaliel_Lopez-Leal*



# QIIME2 using Pair End Mode

All samples were deposited in Bonampack server (132.248.220.35). The work directory is: /space31/PGE/gamlopez/INCAN

We fisrt run Quiime2 in a pair end mode, fisrt we have to activate qiime 

```
source activate qiime2-2018.6
```

Then we have to import the data

```
nohup qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'   --input-path manifest2.csv   --output-path sequences.qza   --source-format PairedEndFastqManifestPhred33 &
```

Here we need two files the manifest file and the mapping file, according to the qiime's manual for pair end mode. The manisfest file looks like:

```
sample-id,absolute-filepath,direction
# Lines starting with '#' are ignored and can be used to create
# "comments" or even "comment out" entries
sample-1,$PWD/some/filepath/sample1_R1.fastq.gz,forward
sample-2,$PWD/some/filepath/sample2_R1.fastq.gz,forward
sample-1,$PWD/some/filepath/sample1_R2.fastq.gz,reverse
sample-2,$PWD/some/filepath/sample2_R2.fastq.gz,reverse
```

You can check the manual here: [https://docs.qiime2.org/2018.11/tutorials/importing/](https://docs.qiime2.org/2018.11/tutorials/importing/)

The mapping file: [https://docs.qiime2.org/2018.2/tutorials/metadata/](https://docs.qiime2.org/2018.2/tutorials/metadata/)

**Note:** In the Pair mode directory (Q2_results) the mapping file is named as metadata.txt



### Demultiplexing the reads

```
qiime dada2 denoise-paired --i-demultiplexed-seqs sequences.qza  --o-table table-dada2_newParameters.qza --o-representative-sequences rep-seqs-dada2.qza  --o-denoising-stats  stats-dada2_newParameters.qza   --p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 200 --p-trunc-len-r 200 --p-n-threads 20
```

To vizualized the repots files:

```
nohup qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv &

nohup qiime feature-table summarize --i-table table-dada2.qza --o-visualization table-dada2.qzv --m-sample-metadata-file metadata.txt &
```



### Build the phylogenetic tree

```
nohup qiime alignment mafft --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs.qza &

nohup qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza &

nohup qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza &

nohup qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza &
```

The outfiles can loaded in iTol to vizualise the tree/trees ([https://itol.embl.de/upload.cgi](https://itol.embl.de/upload.cgi))

### Taxonomic analysis

Before you classified your reads, you have to train the classifier. Here we used GreenGenes data base according to the qiime2 tutorial. Fisrt we download the files to train qiime

```
wget -O "85_otu_taxonomy.txt" "https://data.qiime2.org/2017.7/tutorials/training-feature-classifiers/85_otu_taxonomy.txt"
wget -O "85_otus.fasta" "https://data.qiime2.org/2017.7/tutorials/training-feature-classifiers/85_otus.fasta"
wget -O "rep-seqs.qza" "https://data.qiime2.org/2017.7/tutorials/training-feature-classifiers/rep-seqs.qza"
```

Process the data base and trainning 

```
qiime tools import  --type 'FeatureData[Sequence]'   --input-path 85_otus.fasta --output-path 85_otus.qza

qiime tools import --type 'FeatureData[Taxonomy]'  --source-format HeaderlessTSVTaxonomyFormat --input-path 85_otu_taxonomy.txt --output-path ref-taxonomy.qza

qiime feature-classifier extract-reads  --i-sequences 85_otus.qza --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT --p-trunc-len 100  --o-reads ref-seqs.qza


qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza  --o-classifier classifier.qza
```



Classify your reads and vizualize using qiime2 view ([https://view.qiime2.org/](https://view.qiime2.org/))

```
qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

qiime taxa barplot --i-table table-dada2.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.txt  --o-visualization taxa-bar-plots.qzv
```

### Alpha and beta diversity

```
qiime diversity core-metrics-phylogenetic --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-sampling-depth 19 --m-metadata-file metadata.txt --output-dir core-diversity-phylogenetic_samplingdepth19

qiime emperor plot --i-pcoa bray_curtis_pcoa_results.qza --m-metadata-file /space31/PGE/gamlopez/INCAN/Q2_results/metadata.txt --o-visualization pcoa-visualization.qzv
```

**Note:** using this protocol we lost 90% of the reads in each sample after used dada2 (see dada2 report files). Therefore, we performed the quality filter and merges using Trim_galore, FASTX-Toolkit and PEAR

TrimGalore: [https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

FASTX-Toolkit: [http://hannonlab.cshl.edu/fastx_toolkit/](http://hannonlab.cshl.edu/fastx_toolkit/)

PEAR: [https://cme.h-its.org/exelixis/web/software/pear/](https://cme.h-its.org/exelixis/web/software/pear/)



# QIMME2 usnig Single End Mode

To import the data

```
nohup qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path manifest-SE.csv --output-path single-end-sequence.qza --source-format SingleEndFastqManifestPhred33 &
```

The manifest format looks like:

```
sample-id,absolute-filepath,direction
sample-1,$PWD/some/filepath/sample1_R1.fastq,forward
```

Here we have to denoise the amplicons using dada2:

```
nohup  qiime dada2 denoise-single --i-demultiplexed-seqs single-end-sequence.qza --o-table table-dada2 --o-representative-sequences rep-seqs-dada2 --o-denoising-stats  stats-dada2.qza   --p-n-threads 20 --p-trunc-len 240
```

Check the dada2 table reports files:

```
nohup qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv &

nohup qiime feature-table summarize --i-table table-dada2.qza --o-visualization table-dada2.qzv --m-sample-metadata-file mapping.txt &
```



### Build the phylogenetic tree

```
nohup qiime alignment mafft --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs.qza &

nohup qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza &

nohup qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza &
 
nohup qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza &
```

The outfiles can loaded in iTol to vizualise the tree/trees ([https://itol.embl.de/upload.cgi](https://itol.embl.de/upload.cgi))

### Taxonomic analysis

```
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza  --o-classifier classifier.qza

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

qiime taxa barplot --i-table table-dada2.qza --i-taxonomy taxonomy.qza --m-metadata-file mapping.txt  --o-visualization taxa-bar-plots.qzv
```



### Alpha and beta diversity

```
nohup  qiime diversity core-metrics-phylogenetic --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-sampling-depth 100000 --m-metadata-file mapping.txt --output-dir core-diversity-phylogenetic_samplingdepth100mil &

cd core-diversity-phylogenetic_samplingdepth100mil

qiime emperor plot --i-pcoa bray_curtis_pcoa_results.qza --m-metadata-file /space31/PGE/gamlopez/INCAN/Q2_results-SE/mapping.txt --o-visualization pcoa-visualization.qzv
```



### Alpha rarefaction

```
qiime diversity alpha-rarefaction --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-max-depth 83732 --m-metadata-file mapping.txt  --o-visualization alpha-rarefaction.qzv &
```

**Note:** the value 83732 correspond to the 75% of the reads from the sample with lower reads (after denoised the samples). This was used to include all samples in the test.



