#!/bin/bash
set -u -o pipefail

TOKEN="${TOKEN:-8760468774:AAFDPir-jeXe3oYAxat9eocskvClZ-ebMTw}"
CHAT_ID="${CHAT_ID:-1083936455}"

STATUS="${1:-0}"
if [[ "$STATUS" =~ ^-?[0-9]+$ ]]; then
  shift || true
else
  STATUS=0
fi

if [[ $# -gt 0 ]]; then
  TEXT="$*"
elif [[ "${STATUS}" -eq 0 ]]; then
  TEXT="$(basename "$PWD"): MD completed successfully"
else
  TEXT="$(basename "$PWD"): MD failed (exit ${STATUS})"
fi

curl -s -X POST "https://api.telegram.org/bot${TOKEN}/sendMessage" \
  -d "chat_id=${CHAT_ID}" \
  -d "text=${TEXT}" >/dev/null || true

exit "${STATUS}"
