# Test Reaction SMILES

This file lists diverse reaction SMILES for common transformations. Format: `reactants >> product` (reagents/catalysts omitted).

## Suzuki (C–C)

```
Brc1ccccc1.B(OH)2c1ccccc1>>c1ccc(cc1)c2ccccc2
Brc1ccc(Cl)cc1.B(OH)2c1ccc(F)cc1>>Clc1ccc(cc1)c2ccc(F)cc2
Ic1ccncc1.B(OH)2c1ccccc1>>c1ccc(cc1)c2ccncc2
Brc1ccsc1.B(OH)2c1ccc(OC)cc1>>c1ccc(OC)cc1-c2ccsc2
Brc1ccc(C)cc1.B(OH)2c1cc(C)ccc1>>c1ccc(C)cc1-c2ccc(C)cc2
```

## Buchwald–Hartwig C–N

```
Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Brc1ccc(Cl)cc1.N1CCOCC1>>Clc1ccc(N2CCOCC2)cc1
Brc1ccncc1.NC>>c1cc(NC)ncc1
Brc1ccc(OC)cc1.NCC>>c1ccc(OC)cc1NCC
Brc1cc(C)ccc1.N(C)C>>c1cc(C)ccc1N(C)C
```

## Amide Formation (Acid + Amine)

```
O=C(O)c1ccccc1.NC>>O=C(NC)c1ccccc1
CC(=O)O.NCCC>>CC(=O)NCCC
Clc1ccc(cc1)C(=O)O.N1CCOCC1>>Clc1ccc(cc1)C(=O)N1CCOCC1
O=C(O)C1CCCCC1.Nc2ccccc2>>O=C(Nc2ccccc2)C1CCCCC1
FC(F)(F)c1ccc(cc1)C(=O)O.NC(C)C>>FC(F)(F)c1ccc(cc1)C(=O)NC(C)C
```
