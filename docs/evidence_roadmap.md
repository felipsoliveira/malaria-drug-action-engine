# Evidence-backed scientific roadmap

This document turns literature constraints into implementable, falsifiable extensions of the **Malaria Drug Action Engine (MDAE)**. The goal is not to add mechanisms because they are plausible; each module should have an external observation that can make it fail.

## Priority 1 — Digestive-vacuole exposure v2

### Why

The current chloroquine/PfCRT work already treats digestive-vacuole exposure mechanistically, but a pure weak-base ion-trapping picture is not sufficient.

PfCRT-mediated chloroquine transport is saturable and variant-dependent. Structural work supports a central drug-binding cavity and pH-/membrane-potential-dependent transport. PfMDR1 and PfCRT can also redistribute multiple antimalarials between digestive vacuole and cytosol.

A particularly useful constraint comes from V-ATPase perturbation: severe digestive-vacuole deacidification from roughly pH 4.9 to ~6.0 produced a much smaller reduction in chloroquine accumulation than simple ion trapping predicts. In one experiment, simple partitioning would predict >150-fold loss of accumulation, whereas the observed reduction after pharmacological V-ATPase inhibition was only about threefold after 30 min; chloroquine susceptibility was also largely preserved.

### Proposed model

Extend the exposure layer from

```text
pH partitioning -> local free drug
```

to

```text
passive/pH-dependent partitioning
+ saturable PfCRT efflux
+ PfMDR1-dependent redistribution
+ binding/retention terms
+ variant-specific transport
-> local free drug
```

### Minimum implementation

- explicit Michaelis-Menten PfCRT transport;
- transporter variant object with kinetic parameters and provenance;
- optional PfMDR1 transport term;
- explicit drug-retention term rather than forcing all accumulation into pH partitioning;
- uncertainty on pH, transport parameters and retention;
- tests for limiting cases: no transporter, zero gradient, saturating substrate, and WT/resistant variant swaps.

### Falsifiable benchmark

A future benchmark should require the model to reproduce **both**:

1. saturable/variant-dependent PfCRT transport; and
2. the fact that strong vacuolar deacidification does **not** collapse chloroquine accumulation by the amount predicted by a simple ion-trapping-only model.

### Core literature

- Summers RL et al. *Diverse mutational pathways converge on saturable chloroquine transport via the malaria parasite's chloroquine resistance transporter.* PNAS (2014). https://doi.org/10.1073/pnas.1322965111
- Kim J et al. *Structure and drug resistance of the Plasmodium falciparum transporter PfCRT.* Nature (2019). https://doi.org/10.1038/s41586-019-1795-x
- Shafik SH et al. *Mechanistic basis for multidrug resistance and collateral drug sensitivity conferred to the malaria parasite by polymorphisms in PfMDR1 and PfCRT.* PLoS Biology (2022). https://doi.org/10.1371/journal.pbio.3001616
- Alder A et al. *The role of Plasmodium V-ATPase in vacuolar physiology and antimalarial drug uptake.* PNAS (2023). https://doi.org/10.1073/pnas.2306420120

---

## Priority 2 — Stage-specific artemisinin survival

### Why

A standard long-duration IC50 is a poor phenotype for artemisinin partial resistance. The resistance phenotype is highly stage-dependent.

The RSA0-3h protocol provides an unusually clean benchmark:

- tightly synchronized **0–3 h rings**;
- **700 nM DHA**;
- **6 h exposure**;
- washout and recovery;
- survival measured **66 h after drug removal**;
- survival defined relative to the matched vehicle control.

In the original study, median in-vitro RSA0-3h survival was **0.23%** in fast-clearing isolates and **10.88%** in slow-clearing isolates.

Artemisinin activation is strongly linked to haem. Early rings can draw on parasite haem biosynthesis, whereas later stages rely much more strongly on haemoglobin digestion. K13-associated resistance is linked to reduced haemoglobin endocytosis, which provides a mechanistic bridge from genotype/cell state to reduced drug activation.

### Proposed model

Add parasite age as an explicit state variable:

```text
parasite age
-> haem source / Hb uptake
-> artemisinin activation
-> distributed proteotoxic damage
-> survival probability
```

with an optional K13/endocytosis modifier.

### Minimum implementation

- age bins or continuous hours post-invasion;
- pulse exposure function for DHA;
- activation term driven by stage-dependent haem availability;
- K13/endocytosis modifier;
- survival output at the RSA readout time;
- tests showing that moving the same DHA pulse from 0–3 h to a later age changes predicted survival.

### Falsifiable benchmark

Reproduce the RSA protocol computationally and test whether the model can separate low-survival and high-survival phenotypes without fitting to a conventional 48–72 h IC50.

### Core literature

- Witkowski B et al. *Novel phenotypic assays for the detection of artemisinin-resistant Plasmodium falciparum malaria in Cambodia.* Lancet Infectious Diseases (2013). https://doi.org/10.1016/S1473-3099(13)70252-4
- Wang J et al. *Haem-activated promiscuous targeting of artemisinin in Plasmodium falciparum.* Nature Communications (2015). https://doi.org/10.1038/ncomms10111
- Birnbaum J et al. *A Kelch13-defined endocytosis pathway mediates artemisinin resistance in malaria parasites.* Science (2020). https://doi.org/10.1126/science.aax4735

---

## Priority 3 — Compound property and provenance layer

### Why

Exposure calculations need molecule-specific inputs that are not interchangeable with docking scores. A standardized experimental toolbox for antimalarial compounds reports quantities directly useful to PK/exposure models:

- pKa;
- logD7.4;
- biorelevant solubility;
- permeability;
- plasma/media fraction unbound;
- blood-to-plasma partitioning;
- intrinsic clearance.

The same work also demonstrates that highly lipophilic compounds can produce measurement artifacts through non-specific adsorption, membrane retention, extreme protein binding or apparent metabolic stability.

### Proposed schema

Each property should carry:

```text
value
unit
uncertainty
kind = measured | derived | fitted | assumed
assay / conditions
primary source
notes / limitations
```

### Minimum implementation

Create a `CompoundProperties` structure and keep it independent from target-binding parameters. Exposure modules should consume these properties through explicit interfaces.

### Core literature

- Charman SA et al. *An in vitro toolbox to accelerate anti-malarial drug discovery and development.* Malaria Journal (2020). https://doi.org/10.1186/s12936-019-3075-5

---

## Priority 4 — Apicoplast delayed-death module

### Why

Apicoplast housekeeping inhibitors produce a multi-cycle phenotype that cannot be represented by an instantaneous kill term.

Blood-stage parasites can survive without an apicoplast when supplied with exogenous IPP, demonstrating that isoprenoid precursor production is the critical blood-stage output. Delayed death proceeds through loss of isoprenoid products and protein prenylation, impaired intracellular trafficking, digestive-vacuole dysfunction and growth arrest in the following intraerythrocytic cycle.

Useful experimental constraints include:

- first ~48 h cycle can proceed with little obvious growth defect;
- major lethal phenotype emerges during the second cycle (roughly 48–96 h);
- IPP rescue can maintain parasites despite organelle loss;
- full rescue in the foundational chemical-rescue work used **200 µM IPP**, with partial rescue appearing at lower concentrations;
- loss of prenylation and trafficking defects emerge before terminal arrest.

### Proposed minimal state model

```text
apicoplast function
-> IPP pool
-> FPP/GGPP availability
-> protein prenylation capacity
-> vesicular trafficking / digestive-vacuole function
-> replication competence
```

The model must span multiple parasite cycles.

### Falsifiable benchmark

A valid model should reproduce:

1. delayed rather than immediate death after apicoplast-housekeeping inhibition;
2. rescue by IPP;
3. the ordering of causal events: isoprenoid depletion before prenylation/trafficking failure and growth arrest.

### Core literature

- Yeh E, DeRisi JL. *Chemical Rescue of Malaria Parasites Lacking an Apicoplast Defines Organelle Function in Blood-Stage Plasmodium falciparum.* PLoS Biology (2011). https://doi.org/10.1371/journal.pbio.1001138
- Kennedy K et al. *Delayed death in the malaria parasite Plasmodium falciparum is caused by disruption of prenylation-dependent intracellular trafficking.* PLoS Biology (2019). https://doi.org/10.1371/journal.pbio.3000376

---

## Priority 5 — Follow-up to PfDHFR benchmark 2b

Benchmark 2b established a useful negative result: static GNINA/Vina scores did not recover the experimental pyrimethamine-resistance ordering across the PfDHFR variant panel.

The next computational question should therefore **not** be another re-tuned docking score. A more defensible direction is to test whether explicit conformational sampling and/or alchemical free-energy calculations recover variant-dependent affinity changes.

Before running such a benchmark:

- define the exact variant subset;
- freeze protonation, cofactor, structural preparation and force-field choices;
- pre-register the error metric and acceptance threshold;
- decide how experimental Ki-derived uncertainty is represented;
- include at least one held-out comparison if the data permit.

No PASS threshold should be chosen after seeing the free-energy results.

---

## Implementation order

1. Digestive-vacuole exposure v2.
2. Artemisinin RSA0-3h benchmark.
3. Compound-property schema + uncertainty/provenance.
4. Apicoplast delayed-death module.
5. Alchemical follow-up to PfDHFR benchmark 2b.
6. Integrate cross-drug PfCRT/PfMDR1 collateral sensitivity once the individual transport components are validated.

## Rule for every new module

Every module should enter the repository with:

- a scientific question;
- explicit assumptions and units;
- primary-source provenance;
- a reference or held-out validation target;
- uncertainty where the data support it;
- automated tests;
- a declared failure mode.

A negative result is a valid result.
