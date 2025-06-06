FROM ghcr.io/prefix-dev/pixi:0.40.0 AS build

# copy source code, pixi.toml and pixi.lock to the container
COPY . /app
WORKDIR /app
# run the `install` command (or any other). This will also install the dependencies into `/app/.pixi`
RUN pixi install
# Create the shell-hook bash script to activate the environment
RUN pixi shell-hook -e default > /shell-hook.sh

# extend the shell-hook script to run the command passed to the container
RUN echo 'exec "$@"' >> /shell-hook.sh

FROM ubuntu:24.04 AS production

# Combine apt operations to reduce image layers and cleanup in the same layer
# to minimize size
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    curl \
    sudo \
    zip \
    unzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*


# only copy the production environment into default container
# please note that the "prefix" (path) needs to stay the same as in the build container
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /shell-hook.sh /shell-hook.sh

# Add shell hook sourcing to multiple initialization files to ensure it's loaded
RUN echo 'source /shell-hook.sh' >> /etc/bash.bashrc && \
    echo 'source /shell-hook.sh' >> /root/.bashrc && \
    echo 'source /shell-hook.sh' >> /etc/profile && \
    mkdir -p /etc/profile.d && \
    echo 'source /shell-hook.sh' > /etc/profile.d/pixi_env.sh


WORKDIR /data

# set the entrypoint to the shell-hook script (activate the environment and run the command)
# no more pixi needed in the default container
ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]
