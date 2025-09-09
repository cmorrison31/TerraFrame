#!/usr/bin/env bash
# Usage: ./video2gif.sh input.webm

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 input_video"
    exit 1
fi

INPUT="$1"
BASENAME="${INPUT%.*}"
PALETTE="${BASENAME}_palette.png"
OUTPUT="${BASENAME}.gif"

# Parameters (tweak to taste)
FPS=15
SCALE=640

# Generate palette
ffmpeg -y -i "$INPUT" -vf "fps=$FPS,scale=$SCALE:-1:flags=lanczos,palettegen" "$PALETTE"

# Create gif using palette
ffmpeg -y -i "$INPUT" -i "$PALETTE" -filter_complex "fps=$FPS,scale=$SCALE:-1:flags=lanczos[x];[x][1:v]paletteuse" "$OUTPUT"

echo "Created $OUTPUT"

