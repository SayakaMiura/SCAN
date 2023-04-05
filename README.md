# SCAN_v0.1

(Copyright 2023, Authors and Temple University; see license below)
Updated April 4, 2023

## Description
SCAN improves inferred mutation trees by assessing their reliabilities. The SCAN program was developed by Tenzin Dolker. See Miura et al. (ref. 1) for the detail. You are free to download, modify, and expand this code under a permissive license similar to the BSD 2-Clause License (see below). 

## Dependencies
* python 3 (v3.8.10 was tested)
* python packages: 
 scipy, biopython, numpy, graphviz, pydot, matplotlib
 
 Note: If the installation of these python packages is not easy, you may want to use Anaconda for Python 3 (https://www.anaconda.com/distribution/). Or you can try python3 -m pip install [package name].

## How to use SCAN

`python RunSCAN.py [inferred mutation tree file] [cell order file] [output directory] [inferred cell sequences file] `

#### input 
- inferred mutation tree file (GV format)
 
List (1) the node (mutation ID), (2) the ancestor descendant relationship of mutations, and (3) cell sequence attachement node. Mutation ID is the position of a mutation in the cell sequences file. Cell ID should be named by, ‘s’ number, e.g., s10. The number represent the line of the cell in cell order file. For example,
 
digraph G {
 
node [Mutation ID];

1 -> 2;

1 -> 3;

2 -> 4;

…

1 -> s31;

3 -> s8;

3 -> s83;

…

}
 
- cell order file
 
The line number should correspond to the cell ID in the mutation tree file. For example,
 
Cell1

Cell2

…

- Inferred cell sequences file (FASTA format)
 
"T": Mutant allele

"A": Wild-type allele

"?": Missing base


#### output files
 - Mutation tree with SCAN scores (png and dot file)

- Predicted clone sequences (meg file)
 
- Clone annotation (txt file)
 

#### example
 To run example datasets in Example directory, run the following command.
 
`python RunSCAN.py Example\MutTree.gv Example\CellOrder.txt Test Example\CellSeq.fasta`

 The output files will be stored in "Test" directory.



## Reference:
[1] Sayaka Miura, Tenzin Dolker, Maxwell Sanderford, and Sudhir Kumar. Improving cellular phylogenies through integrated use of mutation order and optimality principles (2023) Submitted to Cancers

## Copyright 2023, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
