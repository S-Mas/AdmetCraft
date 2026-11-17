#!/usr/bin/env bash
set -euo pipefail

export PATH="$HOME/.cargo/bin:$PATH"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

cargo build --workspace --release --offline

mkdir -p "$HOME/.local/bin"
mkdir -p "$HOME/.local/share/applications"
mkdir -p "$HOME/.local/share/icons/hicolor/scalable/apps"

install -m 755 "target/release/admetcraft" "$HOME/.local/bin/admetcraft"
install -m 755 "target/release/admetcraft-cli" "$HOME/.local/bin/admetcraft-cli"
install -m 644 "packaging/admetcraft.desktop" "$HOME/.local/share/applications/admetcraft.desktop"
install -m 644 "packaging/admetcraft.svg" "$HOME/.local/share/icons/hicolor/scalable/apps/admetcraft.svg"

echo "Installed admetcraft and admetcraft-cli into $HOME/.local/bin"
