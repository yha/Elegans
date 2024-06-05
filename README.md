# C. Elegans trajectory and shape analysis.
Includes:
 - Loading of coordinate data and stage timing data
 - Loading and caching of cropped videos
 - Contour computation and storage
 - Midline computation and storage/caching, including head and tail detection

The analysis pipeline was tested with julia version 1.9.3 on Windows, and also with remote workers on a Linux machine
for the contours and midline computation.

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

## Computing contours and midlines
After importing a set of experiments, the contours and midlines are computed by scripts at `scripts/contours_and_midlines`.
Both scripts are configured by editing `scripts/contours_and_midlines/args.toml` Modify the paths listed thereunder `[paths]` 
to point to the experiments, stage data, and desired output paths.
For running the analysis locally, no changes are needed under the `[dist]` section (the default `dist.n_remote = 0` means analysis runs locally. To run analysis remotely, see additional setup below), but you may want to modify the number of worker processes, according to available memory, by setting the value of `n_local`.

After configuring `args.toml`, running the script `scripts/contours_and_midlines/contour-store.jl` will install required packages and
start the contour computation, parallelized among workers. e.g.,
```
julia> include("scripts/contours_and_midlines/contour-store.jl")
[...]
Progress available at localhost:8108 (copied to clipboard)
```
Visit the address listed on the last line to view the progress of workers.
Contours are computed with two different sets of parameters: one for L1, the other for L2-adulthood. Each worker processes one well
with one set of parameters at a time.
Computed contours are stored at the path defined by `paths.contours_dir` at `scripts/contours_and_midlines/args.toml`. There should
be two contour files per well – one for each parameter set.

Midlines computation works similarly, by running the script `scripts/contours_and_midlines/midline-store.jl`. It shares the same configuration file `args.toml`. Computed midpoints are stored at the path defined by `paths.midpoints_dir` in `args.toml`.

### Computing contours and midlines with remote workers

To use remote workers, some setup steps are required on the server and local machine.
The server and local machine should work on the same directories through a network share.
Setup a shared directory on the server accessible on the local machine, where all input file are available. List its path on the remote in `args.toml` under `dist.remote_root`, and its path on the local machine under `paths.local_root`. List the remote's address under `dist.remote`. E.g., if `/home/user/` on `server` is shared and mapped locally to `U:`, your `args.toml` should read
```toml
[paths]
local_root = "U:/"
# ...
[dist]
remote = "server"
remote_root = "/home/user/"
```
Julia should be available on the remote server at `$remote_root/julia-$VERSION/bin/julia`.
Clone or download the repository https://github.com/yha/Elegans/ to the remote. The environment on the remote machine need to be instantiated in advance of running the scripts using `]activate Elegans` followed by `]instantiate`, e.g.:
```
elegans-pipeline> git clone https://github.com/yha/Elegans.git
Cloning into 'Elegans'...
[...]
elegans-pipeline> ls
Elegans
elegans-pipeline> ~/julia-1.9.3/bin/julia
[~] ~/julia-1.9.3/bin/julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.3 (2023-08-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(@v1.9) pkg> activate Elegans/
  Activating project at `~/elegans-pipeline/Elegans`

(Elegans) pkg> instantiate
```
On the local machine, list the location of the repository under `dist.remote_project_root` in `args.toml`, and set `n_remote` to the number of remote workers requires:
```toml
[dist]
# ...
n_remote = 4
# ...
remote_project_root = "path/to/Elegans" # resolved relative to `dist.root`
```

You should also set up password-less ssh access from the local machine to the remote.