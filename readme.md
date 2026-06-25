**Requirements:**

* [BLAT](https://genome.ucsc.edu/goldenPath/help/blatSpec.html)
* python >3.1
* Cython cLev
* Biopython
* R and RStudio
	* dplyr
	* tidyr
	* reshape2
	* magrittr
	* purrr
	* stringdist
	* Biostrings

**1. ./BLAT\_bcs.sh**

Ran BLAT alignments on spkk20, skpp15, and 384 and 1536 barcodes sets (see *seq_files* folder). These were aligned against the orignial [Elledge 25mer 240k](https://elledge.hms.harvard.edu/?page_id=638) barcode set (*seq_files/bc25mer.240k.fasta*). 

**2. python Generate\_barcode\_csv.py**

Nora's code to parse alignments BLAT pretty files (in folder *BCs_psl*) and map all barcodes to original Elledge set.
[https://github.com/PlesaLab/Barcodes\_from\_Elledge](https://github.com/PlesaLab/Barcodes_from_Elledge)

This outputs: *bc25mers\_and\_subsets\_barcodes.csv*

**3. newBCs.R**

This gets all unused Elledge barcodes and outputs them to: *unused\_elledge\_25mers.fasta*

**4. python 1\_12mer\_barcodedesign\_v5\_noLev.py**

This screens for:

output\_primer\_length = 12

cutoff\_GC\_min = 45

cutoff\_GC\_max = 55

lowtemp = 38

hightemp = 44

selfdimercutoff = 4

Restriction sites screened:

* BtsI 
* BspQI
* BsmAI
* BsrDI
* BstNBI
* BsaI
* BbsI
* BsmBI
* BtgZI
* Esp3I
* PaqCI

This saves output:

*filt\_prim\_12nt\_Tm\_38\_44\_GC\_45\_55\_SD\_2.fasta*

with 184,146 sequences.

**5. filt_lev.R**

This is used to filter out sequences closer than Levenshtein distance 3. Note that the data is split into chunks to avoid saturating memory. This operation is VERY memory intense.

This results in 15,500 sequences.

This was then subset into the first 12,288 to create the new bead barcode set.

