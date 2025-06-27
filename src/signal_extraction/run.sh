#! /usr/bin/bash

FMIN=100
FMAX=500

DIR="/media/stagiaireNFS4/maya"
INPUT="$DIR/INTENSO/PADA-1"
OUTPUT="$DIR/results"

echo "Parameters"
echo "  FMIN = $FMIN"
echo "  FMAX = $FMAX"
echo " INPUT = $INPUT"
echo "OUTPUT = $OUTPUT"
echo ""

for folder in "$INPUT"/D1336R*; do
  [ -d "$folder" ] || continue

  ID=$(basename "$folder" | xargs)

  if [[ "$ID" =~ ^D1336R([0][1-9]|[1-7][0-9]|80)$ ]]; then
    echo ">>> Submitting job for sample: $ID"
    echo ">>> BAM will be: $INPUT/$ID/${ID}_dup.bam"

    qsub \
      -q mem48G \
      -l nodes=1:ppn=8 \
      -v ID="$ID",INPUT="$INPUT",OUTPUT="$OUTPUT",FMIN="$FMIN",FMAX="$FMAX",DIR="$DIR" \
      process.sh
  else
    echo "Skipping: $ID (not in R01–R80)"
  fi
done
