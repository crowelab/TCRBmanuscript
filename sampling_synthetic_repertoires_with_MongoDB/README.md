# Sampling Synthetic Repertoires with MongoDB

## Setup

Import a collection of synthetic sequences into MongoDB from a CSV file that includes the headers "v_family,d_family,j_family,cdr3,cdr3_aa_len". It is recommended to use a powerful server to host the synthetic repertoire collection. MongoDB's sampling routines must scan the entire set of matching documents to randomize the selections, so disk speed and the ability to keep most of the collection in memory are derirable.  Creating an index on v_family, d_family, j_family, and cdr3_aa_len fields will dramatically reduce lookup times.

From experimentally derived VDJ+CDR3 length distributions, create a CSV file that contains the VDJ triple separated by underscore (TRBV5-5_TRBD2_TRBJ2-3) followed by a comma-separated list of sample sizes for each CDR3 length (from 0-39 AA residues, inclusive). An example is included for clarity.

## Running

Modify the included MongoDB script to match the desired `num_samples` and the input file (line 35) from the default "example.csv". If your collection is not called `synthetics` the reference to `db.synthetics` on line 3 will also need to be altered.

Then run the following command, replacing `[database name]` with the database where you imported the synthetic sequences. This could be sim1, sim2, or sim3, for example.

`mongo [database name] sample-from-distribution.js > out.dat`

`out.dat` will contain data in the clonotype.dat format used in the other scripts distributed with this GitHub collection. MongoDB imposes a limit of 16MB of memory usage per document, so it was not possible to guarantee a given sampling would fit into that constraint. Therefore, the samples are written below headers "@@synth1.dat", "@@synth2.dat", ... "@@synth1000.dat". Split `out.dat` on these headers and write to individual files.

## Included files

### example.csv
An example of the input for the script that directs it what germline gene families and CDR3 lengths to sample from.

### sample-from-distribution.js
The MongoDB script itself.
