# Installation

This guide walks you through the installation process for Anatalyst.
ALL OF THIS IS WIP AND WILL NEED TO BE UPDATED


## Installation Methods

### Method 1: Using docker (Recommended)

```bash
# Pull the Docker image
docker pull ghcr.io/yourusername/downstream-scrnaseq:latest
```

## Using devcontainer (For Development)

If you're using Visual Studio Code, a devcontainer configuration is provided in the `.devcontainer` folder. This allows you to develop within a container that has all dependencies pre-installed.

1. Install the Remote Development extension pack in VS Code
2. Open the repository folder in VS Code
3. When prompted, click "Reopen in Container"
4. VS Code will build the container and provide you with a fully configured development environment


## Troubleshooting

### Common Issues

#### R Bridge Connection Issues

If the R bridge fails to connect:

```
Error in R bridge: /bin/sh: 1: Rscript: not found
```

Make sure R is installed and the `Rscript` executable is in your PATH.

#### Memory Errors with Large Datasets

If you encounter memory errors when processing large datasets:

```
MemoryError: Unable to allocate array with shape...
```

Try increasing the memory limit in your configuration file:

```yaml
pipeline:
  r_memory_limit_gb: 16  # Increase this value
```

For additional issues, please check the [GitHub Issues](https://github.com/) (dead link for now) page or submit a new issue.

