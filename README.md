# 3D Mixing Model Workflow

MATLAB implementation associated with the paper:

**A 3D inverse design framework for compositional accuracy in DED-processed functionally graded materials**

This repository contains the MATLAB scripts and CSV input files used for the forward dilution model, inverse composition correction, and deviation analysis presented in the paper.

## Requirements

- MATLAB R2024b or later
- MATLAB Optimization Toolbox

## Repository contents

### MATLAB scripts

- `forward_model.m` — computes the output composition from a prescribed input material distribution and toolpath.
- `inverse_model.m` — computes the corrected input composition required to reproduce a prescribed target composition after dilution.
- `deviation_plot.m` — plots material compositions and evaluates deviation from the target or input-correction magnitude.

### Input files

The following material-distribution and toolpath files are included:

```text
input_material_distribution_remelting.csv
toolpath_remelting.csv

input_material_distribution_dragging.csv
toolpath_dragging.csv

input_material_distribution_3d_blisk.csv
toolpath_3d_blisk.csv
```

These correspond to the cases shown in Figs. 2, 3, and 4 of the paper, respectively.

## Usage

### Forward model

Set the input material distribution and toolpath at the top of `forward_model.m`:

```matlab
input_material_file = 'input_material_distribution_3d_blisk.csv';
toolpath_file = 'toolpath_3d_blisk.csv';
```

Run:

```matlab
forward_model
```

The script writes:

```text
output_material_distribution_3d_free.csv
```

This file contains the predicted output composition after dilution.

### Inverse model

Set the target material distribution and toolpath at the top of `inverse_model.m`:

```matlab
material_file = 'input_material_distribution_3d_blisk.csv';
toolpath_file = 'toolpath_3d_blisk.csv';
```

Here, `material_file` is the target composition that the inverse model attempts to reproduce after dilution.

Run:

```matlab
inverse_model
```

The script writes:

```text
required_input_material_distribution_3d_free.csv
validated_output_material_distribution_3d_free.csv
```

- `required_input_material_distribution_3d_free.csv` contains the corrected input composition obtained from the inverse optimization.
- `validated_output_material_distribution_3d_free.csv` contains the predicted output obtained by applying the forward model to the corrected input.

### Deviation and correction plots

Set the target/reference file and the file to be evaluated at the top of `deviation_plot.m`.

For the deviation of the naive output from the target:

```matlab
reference_file = 'input_material_distribution_3d_blisk.csv';
evaluated_file = 'output_material_distribution_3d_free.csv';
plot_mode = "deviation";
```

For the deviation of the corrected output from the target:

```matlab
reference_file = 'input_material_distribution_3d_blisk.csv';
evaluated_file = 'validated_output_material_distribution_3d_free.csv';
plot_mode = "deviation";
```

For the magnitude of the input correction:

```matlab
reference_file = 'input_material_distribution_3d_blisk.csv';
evaluated_file = 'required_input_material_distribution_3d_free.csv';
plot_mode = "correction";
```

Run:

```matlab
deviation_plot
```

## Paper cases

Use the corresponding material-distribution and toolpath files for each case:

| Paper figure | Case | Material distribution | Toolpath |
|---|---|---|---|
| Fig. 2 | Remelting | `input_material_distribution_remelting.csv` | `toolpath_remelting.csv` |
| Fig. 3 | Dragging | `input_material_distribution_dragging.csv` | `toolpath_dragging.csv` |
| Fig. 4 | 3D blisk | `input_material_distribution_3d_blisk.csv` | `toolpath_3d_blisk.csv` |

For each case:

1. Run `forward_model.m`.
2. Use `deviation_plot.m` to evaluate the naive-output deviation.
3. Run `inverse_model.m`.
4. Use `deviation_plot.m` to evaluate the input-correction magnitude.
5. Use `deviation_plot.m` to evaluate the corrected-output deviation.

## Output files

The scripts use the following output filenames:

```text
output_material_distribution_3d_free.csv
required_input_material_distribution_3d_free.csv
validated_output_material_distribution_3d_free.csv
```

Each new run overwrites output CSV files from the previous run. Save or rename the outputs if results from multiple cases need to be retained.

## Citation

If you use this code, please cite the associated paper:

G. Gargiulo, M. Nydegger, Z. C. Cordero, Y. Harduf,  
**A 3D inverse design framework for compositional accuracy in DED-processed functionally graded materials.**
