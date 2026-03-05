# Computational frameworks for de novo design of allosteric biosensors
### Goals of the project
Cells are able to detect and respond to molecules in their surroundings in several ways. One of them is ligand binding to allosteric transcription factors(TFs). These proteins can bind to a specific molecule, triggering a conformational change in the TF that allows it to either promote or repress the activity of a specific gene. This conformational change upon binding is known as allostery.

In vivo, this phenomenom allows cells for example to express antibiotic resistance genes when an antibiotic is present, promote glucose processing pathways when glucose concentration is high... But when these allosteric TFs are used in genetic engineering to regulate the activity of a reporter gene, they effectively act as biosensors for the molecule they bind. This way, the presence in the medium of a molecule of interest can be detected easily, without relying in complex experimental techniques

Apart from employing already existing TFs, protein engineering aims to increase the repertoire of "bio-sensable" molecules by engineering allosteric TF with novel recognition sites  There are many examples of succesful engineered novel biosensors. However, each of them required specific campaigns of mutational scaning or tailored strategies. This master thesis proposes three strategies towards a general solution for this problem, taking advantage of the recent advances in the field of computational de novo design of proteins. The long term goal of this project is developing a computational framework that given an adequate allosteric chasis and a molecule of interest, returns potential functional candidates capable of biosensing the molecule.
### Challenges
This task is obviously challenging. Despite recent advances in active site engineering (with programs like LigandMPNN and RFdiffusion3), finding candidates able to specifically bind a non native molecule still requires a high throughput screening of thousands of candidates. However, the main hurdle lies in preserving the allostery after modifying the allosteric site of the protein. Allosteric mechanisms are usually encoded within structures in a complex, interconnected manner. By achieving binding of a novel molecule, we might disrupt the allosteric behaviour of the system, rendering the protein essentially useless. Therefore, the main challenge here is reshaping the allosteric pocket while preserving the allostericity. We aim to tackle this issue with three strategies of increasing complexity:
- **Strategy 1**: rebuild the protein around a sphere centred around the bound ligand using RFdifussion3 and LigandMPNN followed by a strict filtering process
- **Strategy 2**: rebuild the whole C terminal domain of HTH-like transcription factors, taking advantage of their hinge-like allosteric mechanism
- **Strategy 3**: use a in-house modified RFdiffusion version to generate backbones conditioned with two different DNA binding domain states(bound/unbound) and the presence and absence of a ligand. Through this multi-state protein design, we aim to obtain novel switchable backbones that respond to the ligand.
## Strategies
### Strategy 1: Allosteric binding pocket reconstruction
Found at allosteric_biosensor_design_prediction/strategy_1
In the simplest case scenario, we propose the redesign
### Strategy 1 + Snakemake
Found at allosteric_biosensor_design_prediction/strategy_1_snakemake

### Strategy 2: Rebuild the C-terminal domain of HTH-like transcription factors
Found at allosteric_biosensor_design_prediction/strategy_2

### Strategy 3: De novo generation of allostery through DNA binding domain repositioning
Found at allosteric_biosensor_design_prediction/strategy_3




