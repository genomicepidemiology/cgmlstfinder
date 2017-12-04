############################################################
# Dockerfile to build SeqSero
# # Define app's environment with a Dockerfile so it can be reproduced anywhere:
############################################################

# Set base image to Python Anaconda
FROM continuumio/anaconda

# File Author / Maintainer
MAINTAINER Jose Luis Bellod Cisneros

# Update the repository sources list
RUN apt-get install debian-archive-keyring
RUN apt-key update

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils\
    aufs-tools \
    automake \
    ncbi-blast+ \
    btrfs-tools \
    build-essential \
    ca-certificates \
    curl \
    debian-archive-keyring \
    dpkg-sig \
    emacs \
    git \
    iptables \
    libapparmor-dev \
    libcap-dev \
    libmysqlclient-dev \
    libsqlite3-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    lxc \
    mercurial \
    openssh-server \
    parallel \
    perl \
    pigz \
    r-base \
    reprepro

# Install services from Github (testing repositories)
# TODO Install deployment versions

RUN mkdir /root/.ssh/
RUN ssh-keyscan bitbucket.org >> /root/.ssh/known_hosts
RUN ssh-keyscan cpanmin.us >> /root/.ssh/known_hosts
RUN echo "    IdentityFile ~/.ssh/id_rsa" >> /etc/ssh/ssh_config

RUN chmod 777 -R /tmp && chmod o+t -R /tmp

# CGE Tools
#################################################


# Install dependencies
#################################################

# Install kma #
RUN git clone "https://bitbucket.org/genomicepidemiology/kma.git" /usr/src/kma

# Download services
RUN git clone --recursive  "https://bitbucket.org/genomicepidemiology/cgmlstfinder.git" /usr/src/cgmlstfinder
RUN find -type f -iname "*.py" -print -exec chmod 775 {} \;

# Add service scripts to path
ENV PATH $PATH:/usr/src/cgmlstfinder
ENV PATH $PATH:/opt/conda/bin

# Set convinience aliases
RUN echo "alias edit='emacs'" >> ~/.bashrc
RUN echo "alias ls='ls -h --color=tty'" >> ~/.bashrc
RUN echo "alias ll='ls -lrt'" >> ~/.bashrc
RUN echo "alias l='less'" >> ~/.bashrc
RUN echo "alias du='du -hP --max-depth=1'" >> ~/.bashrc
RUN echo "alias cwd='readlink -f .'" >> ~/.bashrc

# Add PATH to ~/.bashrc -This way you can access programs via ssh
RUN echo "PATH=$PATH" >> ~/.bashrc

# Set default working directory
RUN mkdir /usr/test
WORKDIR /usr/test