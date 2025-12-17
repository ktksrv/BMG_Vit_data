#!/bin/bash

SOURCE="Latest_shear_backup/STZ.py"

if [ ! -f "$SOURCE" ]; then
    echo "ERROR: $SOURCE not found"
    exit 1
fi

for d in Latest_shear_[0-7]*; do
    if [ -d "$d" ] && [ "$d" != "Latest_shear_backup" ]; then
        cp "$SOURCE" "$d/STZ.py"
        echo "Copied STZ.py → $d"
    fi
done

echo "STZ.py synchronized to all Latest_shear_* directories."