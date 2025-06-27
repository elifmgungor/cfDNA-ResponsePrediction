#! /usr/bin/bash
# Activate conda environment
source /media/jnoirelNFS3/miniconda3/etc/profile.d/conda.sh
conda activate /media/jnoirelNFS3/miniconda3/envs/bioinfo

cat <<EOF
Options:
  FMIN = $FMIN
  FMAX = $FMAX
  DIR = $DIR
  INPUT = $INPUT
  OUTPUT = $OUTPUT
  ID = $ID
EOF

mkdir -p "$OUTPUT/$ID"
SCRIPT="/media/stagiaireNFS4/maya/pipeline/score.py"
echo ">>> Running script on BAM file: $INPUT/$ID/${ID}_dup.bam"

python3 "$SCRIPT" \
  --fmin "$FMIN" \
  --fmax "$FMAX" \
  --input "$INPUT" \
  --output "$OUTPUT" \
  "$ID"