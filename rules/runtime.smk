import os


def abs_path(path_value):
    return os.path.abspath(path_value)


def tool_container(tool_name):
    return config["tools"][tool_name]["container"]


def db_mount_dir(tool_name):
    subdir = config["tools"][tool_name].get("db_subdir", tool_name)
    return os.path.join(abs_path(config["execution"]["host_db_root"]), subdir)


def arts_ref_dir(wildcards):
    subdir = config["tools"]["arts"].get("ref_subdir", "arts")
    return os.path.join(abs_path(config["execution"]["host_db_root"]), subdir, arts_reference_set(wildcards))


def arts_root_dir(wildcards):
    subdir = config["tools"]["arts"].get("ref_subdir", "arts")
    return os.path.join(abs_path(config["execution"]["host_db_root"]), subdir)
