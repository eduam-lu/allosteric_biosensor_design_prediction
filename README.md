# Computational frameworks for de novo design of allosteric biosensors
### Goals of the project
Cells are able to detect and respond to molecules in their surroundings in several ways. One of them is ligand binding to allosteric transcription factors(TFs). These proteins can bind to a specific molecule, triggering a conformational change in the TF that allows it to either promote or repress the activity of a specific gene. This conformational change upon binding is known as allostery. (Figure 1)

In vivo, this phenomenom allows cells for example to express antibiotic resistance genes when an antibiotic is present, promote glucose processing pathways when glucose concentration is high... But when these allosteric TFs are used in genetic engineering to regulate the activity of a reporter gene, they effectively act as biosensors for the molecule they bind. This way, the presence in the medium of a molecule of interest can be detected easily, without relying in complex experimental techniques. (Figure 2)

Apart from employing already existing TFs, protein engineering aims to increase the repertoire of "bio-sensable" molecules by engineering allosteric TF with novel recognition sites  There are many examples of succesful engineered novel biosensors. However, each of them required specific campaigns of mutational scaning or tailored strategies. This master thesis proposes three strategies towards a general solution for this problem, taking advantage of the recent advances in the field of computational de novo design of proteins. The long term goal of this project is developing a computational framework that given an adequate allosteric chasis and a molecule of interest, returns potential functional candidates capable of biosensing the molecule. 
### Challenges
This task is obviously challenging. Despite recent advances in active site engineering (with programs like LigandMPNN and RFdiffusion3), finding candidates able to specifically bind a non native molecule still requires a high throughput screening of thousands of candidates. However, the main hurdle lies in preserving the allostery after modifying the allosteric site of the protein. Allosteric mechanisms are usually encoded within structures in a complex, interconnected manner. By achieving binding of a novel molecule, we might disrupt the allosteric behaviour of the system, rendering the protein essentially useless. Therefore, the main challenge here is reshaping the allosteric pocket while preserving the allostericity. We aim to tackle this issue with three strategies of increasing complexity (Figure 3):
- **Strategy 1**: rebuild the protein around a sphere centred around the bound ligand using RFdifussion3 and LigandMPNN followed by a strict filtering process
- **Strategy 2**: rebuild the whole C terminal domain of HTH-like transcription factors, taking advantage of their hinge-like allosteric mechanism
- **Strategy 3**: use a in-house modified RFdiffusion version to generate backbones conditioned with two different DNA binding domain states(bound/unbound) and the presence and absence of a ligand. Through this multi-state protein design, we aim to obtain novel switchable backbones that respond to the ligand.

<img src="img/Figure 3.png" alt="Figure 3. Strategy showcase">

## Strategies
### Strategy 1: Allosteric binding pocket reconstruction
Found at [Strategy 1 folder](strategy_1/)

In the simplest case scenario, we propose placing the molecule of interest in the allosteric site of the target protein, and rebuild all the residues found within a sphere of radius X armstrongs with RFdiffusion3, following a sequence recovery step with Ligand MPNN and a high throughput filtering of the candidates (Figure 4).

Through this approach, and by testing different radii and parameters, we expect to find a configuration that can produce candidates whose reconstruction was good enough to fit the new ligand and at the same time, respected the allosteric mechanism.

There are several flaws to this approach however. First, given that allosteric mechanisms are deeply encoded within structures, even small changes in the allosteric site can lead to their complete eradication. Furthermore, some allosteric mechanisms are triggered by the specific chemical properties of the ligand. So even if we managed to fit a new ligand and preserve all the elements of allostery, if the new molecule was not as charged or electronegative as the original, the required structural reorganisation is not possible
### Strategy 1 + Snakemake
Found at [Strategy 1 + snakemake folder](strategy_1_snakemake/)

RFdiffusion and structure prediction programs are heavily dependent on GPU computing. The code for strategy 1 is designed to work with 1 GPU only, making it linear, slow and inneficient. For that reason, we adapted the code to work as a snakemake pipeline.

Snakemake is a workflow managing system designed to make reproducible analyses and handle the available resources in an efficient manner. With its rule based syntax, it predetermines the number of jobs that will be required for each step of the pipeline and the resources they require and manages them throughout the whole execution.

Additionally, in this step we also divided RF diffusion and sequence recovery together with the posterior processing into sepparate modules. In the following strategies, the RFdiffusion approach will be completely replaced, while the sequence recovery and candidate selection module will be updated with novel functionalities.
### Strategy 2: Rebuild the C-terminal domain of HTH-like transcription factors
Found at [Strategy 2 folder](strategy_2/)

The second strategy was then designed to address the flaws of the first one and was built upon the snakemake code for strategy 1. We need allosteric mechanisms that A) are not dependent on the ligand properties and B) are resilient even upon modifying the allosteric site. The HTH (Helix-Turn-Helix) family of TFs and their "hinge-like" allosteric mechanism fully meet those conditions. Allosteric HTH TFs posess independent C and N terminal domains that can move freely and relatively to each other in a hinge-like motion. The allosteric site is found at the intersection of these two domains. Binding of the ligand sterically impedes the hinge like motion, fixating the protein in a state that is either compatible with DNA binding (promoting gene expression) or incompatible with DNA binding (repressing gene expression)(Figure 5).

In these systems, allostery is not necessarily hardly encoded in the structure, and the mere presence of the ligand is enough to block the hinge-like motion, regardless of its properties in principle. There are other advantages: because of the independence of the two domains, one of them can be rebuilt around the ligand to generate a proper pocket, leaving the other almost untouched. This gives RFdiffusion a lot more flexibility than in the previous strategy, conditions under which is known to perform better. Also, given the variability found within this family, several chasises can be tested to accomodate ligands of different sizes and properties.

In this strategy, we propose a full redesign around a molecule of interest of C terminal domains of HTH allosteric TFs, followed by sequence recovery by LigandMPNN. Succesful candidates will be selected based on their ability to bind the new ligand but also for their potential to show a hinge-like motion. For prompting this behaviour, we include multi-state design elements in the process, including a multistate version of LigandMPNN that returns sequences accounting for two potential backbones, template conditioning for structure prediction methods and measuring the energy difference between the two potential states (Figure 6).

Apart from the traditional hurdles of de novo protein design, two challenges emerge from this approach:
- HTH TFs functional version is a homodimer, sometimes even a homotetramer. And although the C termini is not directly involved in allostery, it is involved in dimerization. While reconstructing the C terminal domain, we also need to recover a proper dimerization interface, we account for this in our design process.
- The DNA binding region structure has been rarely determined crystalographically and is usually not present in PDB structures for these TFs. Although distant for the reconstruction site, this might be a source of error.

### Strategy 3: De novo generation of allostery through DNA binding domain repositioning
Found at [Strategy 3 folder](strategy_3/)

Finally, we propose an alternative, fully de novo, strategy for developing novel biosensors. Rather than relying on natural backbones, we propose using an in-house version of RF diffusion (under development by Samuel Sellberg), capable of being conditioned with two protein states, to generate multistate compatible backbones. 

Denoising diffusion models essentially transform Gaussian noise into a denoised system, in this case, a viable protein structure. This denoising process is iterative and can be conditioned by context, for example, by providing a protein region that needs to be included in the final design. Through denoising, the cloud of points adapts to the context to accept it, it gathers information from it. We propose that by exposing the noise to two different contexts throughout denoising, we can generate designs compatible with both. This involves sequentially moving the point cloud between contexts after applying x denoising steps, until the noise collapses in a final structure (Figure 7).

For this project, the contexts consist on a natural DNA binding domain in a DNA binding compatible position and in an incompatible position. Also, the ligand should be present in either of the states, so its presence and absence can also condition the design. From this set up, we obtain to obtain designs with features compatible with both states in a manner regulated by the presence and abscence of the ligand.In other words, we expect to design multi-state backbones where allostery might arise (Figure 8).

This is is a high-risk high-reward approach. This backbone conditioned or conscious redesign is newly developed and would be benchmarked through this project. Also, multistate protein design is levels of magnitude more complicated than traditional protein design. But if successful, this method could be extremely powerful and generalizable. Not only one could fit virtually any molecule of interest, one could also choose the promoter system it regulates by selecting an adequate DNA binding domain.



