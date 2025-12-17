#!/bin/bash
set -e

TEMPLATE="Latest_shear_backup/LAMMPS_slurmscript.sh"

if [ ! -f "$TEMPLATE" ]; then
    echo "ERROR: Template not found: $TEMPLATE"
    exit 1
fi

for d in Latest_shear_[0-9]*; do
    [ "$d" = "Latest_shear_backup" ] && continue

    i=${d##*_}
    OUT="$d/LAMMPS_slurmscript.sh"

    sed "s/{{I}}/$i/g" "$TEMPLATE" > "$OUT"

    echo "Updated $OUT (i=$i)"
done

echo "All SLURM files regenerated from template."