# Campaign specifications

The exact inputs with which the design campaigns were launched on the compute
machine, kept here so that a campaign can be repeated without depending on
that machine's file system.

## rbd/

Nanobody campaign against the SARS-CoV-2 RBD.

- `7BEP_chainE_boltzgen_target.cif`: the target exactly as BoltzGen ingests it.
  Chain `Exp`, 197 residues, numbering 1..197.
- `design_spec.yaml`: the design specification. The epitope window
  `binding: 148..163` refers to that 1..197 numbering and corresponds to
  residues 479..494 of the Spike protein (`PCNGVEGFNCYFPLQS`), the core of the
  RBM that includes E484, F486, Y489 and Q493.

Launch command:

```
boltzgen run --protocol nanobody-anything --output run1 --devices 4 \
             --num_designs 10000 --diffusion_batch_size 1 --reuse \
             --budget 100 design_spec.yaml
```

### Relationship with `validation_panel/reference_targets/`

`validation_panel/reference_targets/7BEP_chainE_boltzgen_target.cif` contains
**the same structure** (197 residues, 1556 atoms, maximum coordinate
deviation 0.0000 Å) but with the chain named `E` and Spike's original
numbering (332..528). That variant serves to compare against
PDB structures; the one in `campaign_specs/rbd/` is the one that
`design_spec.yaml` consumes. They are not interchangeable: with the panel variant, the
`chain: id: Exp` does not resolve and the `148..163` window would select
different residues.

### Provenance of the scaffolds

The four nanobody scaffold specifications (`7eow.yaml`,
`7xl0.yaml`, `8coh.yaml`, `8z8v.yaml`) that `design_spec.yaml` references **are
not artifacts of this project**: they are example files from BoltzGen itself,
unmodified. They live in `boltzgen/example/nanobody_scaffolds/`, which this repo
ignores on purpose (see `.gitignore`) because it is a clone of a third-party repository.

They are reproduced by cloning the fork at commit `7ae564b`. Verified: the
local copies are byte-for-byte identical to the ones the campaign consumed on the
compute machine, and `git status` over that directory is clean with respect to the fork.
