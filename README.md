# IgTreeZ
A toolkit for immunoglobulin gene lineage tree analysis.
Hadas Neuman and Ramit Mehr

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

To run IgTreeZ, you will need these python packages:
KKKKDFFDFGD

To use the plot option, you will need this package

To use draw sub-program, you will need the this program
dot version 2.30.1 or higher

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/neumanh/IgTreeZ.git
   ```
2. Download tar.gz file for linux
   ```sh
   npm install
   ```
3. Download zip file for windows
   ```sh
   npm install
   ```

<!-- USAGE EXAMPLES -->
## Usage

IgTreeZ program includes 5 sub-programs:
* mutations
* poptree
* mtree
* filter
* draw   

You can list them using
   ```sh
   igtreez.py -h
   ```
### mutations

*IgTreeZ mutations* counts and profiles the mutations in a repertoire, based on tree topology. The input for this program must include the trees and sequences, which can be given as a Fasta file or as AIRR/Change-O database, or as an AIRR JSON rearrangement scheme, which includes the trees and sequences.

#### Using AIRR/Change-O database

```sh
igtreez.py mutations -n example_name -t ../examples/*nw -d ../examples/F1-control_germ-pass.tab
```

Notes:
The `-t` parameter can get files of folders.
The input database can be in Change-O of in AIRR format. Use the `-dbf` parameter to specifiy the input format using `-dbf airr` or `-dbf changeo`. The default is the (old) Change-O format.
The program uses the database fields (for Change-O / AIRR format):
* SEQUENCE_IMGT / sequence_alignment 
* GERMLINE_IMGT_D_MASK / germline_alignment 
* CLONE /clone_id 
* SEQUENCE_ID / sequence_id 
You can specify different fields names using the parameters:
`-sf`, `-gf`, `cf`, `if`.

If the CDR3_IMGT / cdr3_imgt column exists in the database, and the '--nocdr3' parameter is avoided - the program defines mutations in the CDR3 region as well

The program profiles the mutation region based on [IMGT region definition](http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html).

#### Using Fasta file example
   ```sh
  igtreez.py mutations -n example_name -t ../examples/1004.nw -f ../examples/1004_aligned.fasta
  ```
Notes:
Eche Fasta file represents the sequences of one tree. The number of the trees must be identical to the number of Fasta files. 
Each Fasta file must include sequences in the names of the corresponded tree's nodes, and the germline sequence.
The program assumes the germline sequence is named **naive**. You can define another name using the `-gl` parameter. All Fasta files must include one germline sequence.
All the sequences in one Fasta file must be aligned (= have same length).

#### Using AIRR scheme
```sh
igtreez.py mutations -n example_name -j ../examples/full_schema_dataset_example.json
```

#### Outputs
The IgTreeZ-mutation program creates a CSV files named *IgTreeZ_output/example_name/example_name_mutations.csv*  which include 82 mutation properties. Each line in the file refers to one tree and each column to a different mutation property. The properties are described [here](link)

##### Plots
In addition to the mutations count, the program can generate discriptive plots using the `--plot` parameter. For example:

All the plots are created in the local directory *IgTreeZ_output/example_name/mutation_plots*. The program generates more than 75 plots for each analysis. A short explanation on each plot can be found [here](link)

An example plot that describes....


<!-- LICENSE -->
## License

Distributed under the AGPL3 License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

Hadas Neuman hadas.doron@gmail.com
Prof. Ramit Mehr ramit.mehr@biu.ac.il 


Project Link: [https://github.com/neumanh/IgTreeZ](https://github.com/neumanh/IgTreeZ)


   

