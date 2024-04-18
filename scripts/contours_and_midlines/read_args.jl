using TOML

"""
Parse the `args.toml` file that is used by the contours and midlines scripts,
and resolve paths.
"""
function read_args(local_project_root, args_path = "$local_project_root/args.toml")
    args = TOML.parsefile(args_path)
    paths, dist = args["paths"], args["dist"]

    stage_path   = normpath(joinpath(local_project_root, paths["stage_path"]))

    local_root = paths["local_root"]
    logdir = paths["logdir"]

    remote_root = dist["remote_root"]
    remote_project_root = "$remote_root/$(dist["remote_project_root"])"

    ex_dir         = paths["ex_dir"]
    contours_dir   = paths["contours_dir"]
    midpoints_dir  = paths["midpoints_dir"]

    local_ex_path       = normpath(joinpath(local_root, ex_dir))
    local_contours_path = normpath(joinpath(local_root, contours_dir))

    assert_isdir(path)  = isdir(path)  || error("Not found or not a directory: $path")
    assert_isfile(path) = isfile(path) || error("Not found or not a regular file: $path")

    assert_isdir(local_root)
    assert_isdir(local_project_root)
    assert_isdir(local_ex_path)
    assert_isdir(local_contours_path)
    assert_isfile(stage_path)

    remote = dist["remote"]
    n_local, n_remote = dist["n_local"], dist["n_remote"]


    skip_existing = args["params"]["skip_existing"]
    
    (;
        stage_path, local_root, remote_root, remote_project_root, logdir,
        ex_dir, contours_dir, midpoints_dir, local_ex_path, local_contours_path,
        remote, n_local, n_remote,
        skip_existing
    )
end