#!/bin/bash
# Runs one DZA01.py cycle. Invoked by launchd on a schedule -- see
# deployment/com.rvx.dza01-sonify.plist and the README's "Deployment
# (macOS / unattended runs)" section. Edit the action list/flags below to
# change what each scheduled run does (see README's Command-line reference).
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate
python3 DZA01.py fetch plot sonify
