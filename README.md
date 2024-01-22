# Molecular complexity

Calculate molecular complexity from a V2000 Molfile.

Original column:\
M. von Korff, T. Sander, *Chimia*, **2023**, 77, 258, DOI: [10.2533/chimia.2023.258](doi.org/10.2533/chimia.2023.258).

Original Java implementation:\
[github.com/Actelion/openchemlib](https://github.com/Actelion/openchemlib)

**Note:** Results may vary between the programs. This is most likely due to how chiral molecules are handled. May be fixed in the future? 

## Requirements
- python (3.9.16 tested)
- numpy (1.22.4 tested)
- networkx (3.0 tested)

## Usage 

```
$ python complexity.py mols/alanine.mol
1.6309297535714573
```