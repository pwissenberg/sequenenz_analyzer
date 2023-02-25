# sequenenz_analyzer
## Introduction
This project is university project.
In this project we were provided with a data set. The data set contained several protein prediction-embeddings. In addition, we had for each protein the information if it the protein has an exophysiologial(exo.) function or an endophysiological(endo.) function. 

After conducting a normal EDA, we were trying to analyse the data in a proper way. Therefore, I came up with the idea to analyse all embeddings which are to a prediction for each amino acids. To find then a combination or a specific amino acids which has a great evidence for each amino acid to be exo or endo. This idea tries to find a new way to look at each specific amino acids and the connected properties of the amino acid.

## Wrangling with the embedddings
In addition, I only looked at the following embeddings:
- Amino Acid sequence
- Dssp3 Prediction for each amino acid
- Dssp8 Prediction for each amino acid
- Conservation score for each amino acid
- Membran_tmbed score for each amino acid
- Metall Binding per amino acid
- Nucleic Acid Binding per amino acid
- Small Molecules Binding per amino acid

```
position:                   0   1   2   3   4   5   6   7   8   9   ...
AA:                         V   G   V   L   L   D   I   L   Q   R   ...
dssp3:                      E   L   L   L   H   H   H   H   E   L   ...
dssp8:                      H   H   H   T   T   C   E   E   E   E   ...
conserv_score:              0   3   3   3   4   0   5   0   5   5   ...
membran_tmbed:              S   S   S   o   o   o   o   o   o   o   ...
metal_binding:              M   M   -   M   -   -   -   -   -   -   ...
nucleic_binding:            N   -   -   -   N   -   -   -   -   -   ... 
small_molecules_binding:    -   -   S   -   -   S   -   S   -   S   ...
```

In the next phase, I was constructing for each amino acid a `word` by combining all the the letters under each position. Each `word` describes the amino acid more precisely. With the example above I was constructing the following words for position 0, 1 and 2:
```
0:  VEH0SMN-    |   Translated: Valine,Sheet,alpha-helix,cons:8,Transmembrane-alpha-helix,MetallBinding,NucleicBindingSmallMoleculeBinding
1:  GLH3SM--    |   Tranalated:...
2:  VLH3S--S    |   Translated:...
```
So, I was creating for each of the proteins their set of `amino-acid-word`. Therefore, I was collecting in addition the positions where the words occure and how many the word occured per data set.
## Interpreting the data parts
