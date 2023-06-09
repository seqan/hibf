#!/usr/bin/env bash
set -Eeuo pipefail

sudo sed -i 's@azure.archive.ubuntu.com@mirror.enzu.com@' /etc/apt/sources.list
sudo add-apt-repository --no-update --yes ppa:ubuntu-toolchain-r/ppa
sudo add-apt-repository --no-update --yes ppa:ubuntu-toolchain-r/test
sudo wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | sudo tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc
sudo add-apt-repository --no-update --yes "deb http://apt.llvm.org/jammy/ llvm-toolchain-jammy main"
sudo apt-get update
