#!/usr/bin/env bash
set -Eeuo pipefail

sudo sed -i 's@azure.archive.ubuntu.com@mirror.enzu.com@' /etc/apt/sources.list
sudo add-apt-repository --no-update --yes ppa:ubuntu-toolchain-r/ppa
sudo add-apt-repository --no-update --yes ppa:ubuntu-toolchain-r/test
sudo apt-get update
