# Rice Microbiome (2018)

Preface:  There will be an upload of a master script that contains the entire data processing, ASV assignment, etc. that's much more legible/user friendly.  It's only not uploaded because of a recurring RAM issue.

These are a series of R scripts that were used to process and run analyses on data that was output from a Illumina MiSeq 300bp Paired End Sequencer run.  The usage of R began after the amplicon sequencing variants (ASVs, similar to OTUs, but assigned by DADA2) were assigned. 

These sets of codes were written during my earliest stages of coding and are fairly trashy and lacking (sorry!) - writing code for others and loops weren't understood at the time.  As such, there are some pieces that lack continuity and many places of excess/unnecessary code due to data exploration/experimentation.

### Prerequisites

R and associated packages (microbiome, DADA2, vegan, and more)

### Installation 

Extract files to a directory.  

## Description of Scripts

DADAPrep:         Uploading of ASV(OTU) files, rarefication, creation of various subsets of the data based on abundance and metadata, selection of a core microbiome

Diversity:        Calculation of Shannon Diversity

PCA:              Principle Component Analysis and figures

VennDiagram:      Generation of 4 part Venn Diagram displaying number of unique ASVs based on metadata categories

ClusterAssign:    Assigns ASVs into 1 of 6 clusters based on abundance patterns across metadata categories

CompositionChart: Generates composition charts for both read and ASV abundance (very ugly code)
