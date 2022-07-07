# Containerizing a repository for use with VSCode, Codespaces, or HPC via Singularity

## Add this as a submodule to the project repository

You can make a repository containerized by cloning this repo as a submodule of your project folder. In the project root, 

```
git submodule add https://github.com/drejom/devcontainer_RbioC.git .devcontainer
git commit -am 'Add devcontainer_RbioC module'
```

### Local R packages
If there's a local package cache to use, set the local folder in `.devcontainers.env`:

`R_USR_LIBS=/local/path` 

## Additional features

* Fully set up for R in VSCode; dependencies are included in the container build

* Supports use of a local package cache if available. If not, packages are installed in the container (they'll be lost if not committed to the image)

* Includes some convenience packages such as `fnmate` and `datapasta` out of the box

## Customizing the shell and R

Edit `.zshrc`, `.lintr` and `.Rprofile` in this folder to customise the *Zsh* and *R* environments.

## HPC Support

Slurm (17.11.13); needs munge.key and slurm.conf mounted from the host. 
Use `future` or `clustermq` for parallel computing on multiprocessor or cluster architectures

Works with Singularity; not yet with Docker

## RStudio support

Run the container and point your browser to http://[container_host]:8787 with your username and the password `rstudio`

## Singularity

Run this on the HPC via a custom [ssh command for VSCode](https://github.com/microsoft/vscode-remote-release/issues/3066); use templates for `batchtools` or `clustermq` for use in R (eg, in `targets` pipelines)

See ~/clustermq_[ssh|slurm].sh for submission scripts