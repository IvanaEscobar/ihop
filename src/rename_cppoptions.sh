#!/usr/bin/env bash

set -euo pipefail

# Build sed script
SED_SCRIPT=$(mktemp)

cat << 'EOF' > "$SED_SCRIPT"
s/\belli_mod/belli_mod/g
s/\belliDoNcOutput/belliDoNcOutput/g
EOF
#cat << 'EOF' > "$SED_SCRIPT"
#s/\BELLI_DEBUG/BELLI_DEBUG/g
#s/\BELLI_THREED/BELLI_THREED/g
#s/\BELLI_3D_STATE/BELLI_3D_STATE/g
#s/\BELLI_2D_STATE/BELLI_2D_STATE/g
#s/\BELLI_TENDENCY/BELLI_TENDENCY/g
#s/\BELLI_MULTIPLE_TIMES/BELLI_MULTIPLE_TIMES/g
#s/\BELLI_MULTIPLE_SOURCES/BELLI_MULTIPLE_SOURCES/g
#s/\BELLI_MULTIPLE_RECEIVER_DEPTHS/BELLI_MULTIPLE_RECEIVER_DEPTHS/g
#s/\BELLI_MULTIPLE_RECEIVER_RANGES/BELLI_MULTIPLE_RECEIVER_RANGES/g
#s/\TEST_BELLI_COST/TEST_BELLI_COST/g
#s/\TEST_BELLI_COST_INLOOP/TEST_BELLI_COST_INLOOP/g
#EOF

#cat << 'EOF' > "$SED_SCRIPT"
#s/IHOP_SIZE\.h/BELLI_SIZE.h/g
#s/IHOP_OPTIONS\.h/BELLI_OPTIONS.h/g
#s/IHOP_COST\.h/BELLI_COST.h/g
#s/IHOP\.h/BELLI.h/g
#s/\ALLOW_BELLI/ALLOW_BELLI/g
#s/\BELLI_WRITE_OUT/BELLI_WRITE_OUT/g
#EOF

changed=0

for f in *; do
  [ -f "$f" ] || continue

  sed -i.bak -f "$SED_SCRIPT" "$f"

  if ! cmp -s "$f" "$f.bak"; then
    echo "Updated: $f"
    changed=$((changed+1))
  fi

  rm -f "$f.bak"
done

rm -f "$SED_SCRIPT"

echo "Done. Files changed: $changed"
