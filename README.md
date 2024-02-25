# C. Elegans trajectory and shape analysis.
Includes:
 - Loading of coordinate data and stage timing data
 - Loading and caching of cropped videos
 - Contour computation and storage
 - Midpoints computation and storage/caching, including head and tail detection

## Importing experiments
To import experiment data, after running the MATLAB analysis scripts, several Julia scripts need to be run
on (copies of) the MATLAB analysis output directories.
These scripts are currently designed to be run inside VSCode, by manually editing input argument assignments 
(and possibly some code) at the top of the script, similarly to the MATLAB scripts. 
They are also divided into cells (separated by `##`), so they can be run one cell at a time if desired to examine 
intermediate results.

The import scripts are at `scripts/experiment_import`, and have their own Julia environment which should be activated 
first in a fresh julia REPL.
After installing vscode and the vscode julia extension, start a Julia REPL inside VSCode (`Ctrl+Shift+p` for the command palette
then choose "Julia: Start REPL", or `Alt+j,Alt+o`).
To activate the environment you can run `]activate scripts/experiment_import` at the REPL running in the `Elegans` project directory 
(hitting `]` in the REPL enters "`pkg` mode", where you type the rest of the command). Alternatively, in VSCode you can right-click 
the directory and select "Julia: activate this environment" (this option is only available when a Julia REPL has been started).
When doing this for the first time on given machine, the required packages are probably not installed yet, so you need to also run
`]instantiate`.

The scripts should generally be run in the following order:

- `denest_wells.jl`: 

    MATLAB coordinate analysis output from each experiment is sometimes divided into subdirs when it was run in several batches.
    This scripts flatten the directory structure to so each well is at `<experiment>/<well>`.

    Modify the assigment to `root` at the top of the script. 
    Comment out the line at the bottom to do a real run rather than just listing which files are moved.

- `save_to_single_file.jl`

    Saves coordinates, video file frame-number ranges, and size information to a single file per well, named `coords_and_size.jld2`.

    Modify the assigment to `root` and the regular expression for experiment names (or the whole condition) in the following line.
    Below that, you may list wells to skip in the `badwells` variable.

- `import_stages.jl`

    This script operates on the "analyzed" experiment data (the output of the three MATLAB scripts starting with `separate_to_developmental_stages.m`), as well as coordinate and video data.
    It imports stage data (which frame each stage starts with for each well) from MATLAB output into the `stages.toml` file stored inside the `Elegans` repository.

    Modify `coord_data_root` to point to coordinate and video data and `analyzed_root` to analyzed data, and the `expname` function logic according to how experiment directories are named in this analysis (this function should extract experiment name from dirname).

