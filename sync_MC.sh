#!/bin/bash
set -e

TEMPLATE="Latest_shear_backup/MC_test.py"

if [ ! -f "$TEMPLATE" ]; then
    echo "ERROR: Template not found: $TEMPLATE"
    exit 1
fi

for d in Latest_shear_[0-9]*; do
    [ "$d" = "Latest_shear_backup" ] && continue

    i=${d##*_}
    OUT="$d/MC_test_$i.py"

    sed "s/{{I}}/$i/g" "$TEMPLATE" > "$OUT"

    echo "Generated $OUT (current_proc=$i)"
done

echo "All MC_test files synchronized."