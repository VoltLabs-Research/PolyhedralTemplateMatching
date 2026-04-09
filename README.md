# PolyhedralTemplateMatching

`PolyhedralTemplateMatching` classifies atoms using PTM and exports the reconstructed state consumed by downstream DXA-compatible tools.

## CLI

Usage:

```bash
polyhedral-template-matching <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--crystalStructure <type>` | No | Input crystal structure: `SC`, `FCC`, `HCP`, `BCC`, `CUBIC_DIAMOND`, `HEX_DIAMOND`. | `FCC` |
| `--rmsd <float>` | No | RMSD threshold for PTM. | `0.1` |
| `--dissolveSmallClusters` | No | Mark small clusters as `OTHER` after clustering. | `false` |
| `--help` | No | Print CLI help. | |

## Build With CoreToolkit

```bash
cd /path/to/voltlabs-ecosystem/tools/CoreToolkit
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/StructureIdentification
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/CommonNeighborAnalysis
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/PolyhedralTemplateMatching
conan create . -nr
```
