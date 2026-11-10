# Docker image to provide a Ubuntu-24.04 environment with debugging tools for ECOGEN: gdb, valgrind, addressSanitizer, tmpi,...
#   - to load the prebuilt image:
#       docker run -it 7g8efcehpff/ecogen_dbg
#   - to build a new multi-platform image from this dockerfile:
#       docker buildx build --platform linux/amd64,linux/arm64 -t <image_name> .
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        binutils \
        build-essential \
        bzip2 \
        ca-certificates \
        clang \
        clang-tidy \
        clang-tools \
        cmake \
        cmake-curses-gui \
        curl \
        dirmngr \
        dpkg-dev \
        fakeroot \
        file \
        g++ \
        gcc \
        gdb \
        gfortran \
        git \
        gnupg \
        less \
        make \
        libopenmpi-dev \
        openssh-client \
        pkg-config \
        python3 \
        python3-pygments \
        python3-yaml \
        reptyr \
        software-properties-common \
        tmux \
        valgrind \
        vim \
        wget \
        xz-utils \
    && rm -rf /var/lib/apt/lists/*

RUN curl https://raw.githubusercontent.com/Azrael3000/tmpi/master/tmpi -o /usr/local/bin/tmpi && chmod +x /usr/local/bin/tmpi

RUN <<'EOF'
cat > /root/.tmux.conf <<'TMUX'
set -g default-terminal "screen-256color"

set -g history-limit 5000
set -g mouse on

bind K kill-session

# easily toggle synchronization (mnemonic: e is for echo)
# sends input to all panes in a given window.
bind e setw synchronize-panes on
bind E setw synchronize-panes off

TMUX
EOF
