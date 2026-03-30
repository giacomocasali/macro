#!/bin/bash

# Official code language: English [cite: 2026-02-26]

echo "-------------------------------------------------------"
echo "ROOT Analysis Directory Cleaner"
echo "-------------------------------------------------------"

# Always run from the directory where this script lives,
# regardless of where it is called from.
cd "$(dirname "$0")"

# 1. Automatic removal of ACLiC generated files and shared libraries
echo "Step 1: Removing compiled libraries and temporary files..."
rm -vf *.so *.d *.pcm *_ACLiC_dict* *_rdict.pcm

# 2. Automatic removal of Windows Zone.Identifier files
echo "Step 2: Removing :Zone.Identifier files..."
find . -name "*:Zone.Identifier" -type f -delete
echo "Done."

# 3. Interactive removal of PNG plots
echo "-------------------------------------------------------"
read -p "Do you want to delete all .png files in this directory? (y/n): " choice

case "$choice" in
  y|Y )
    echo "Deleting all PNG files..."
    rm -vf *.png
    echo "Cleanup complete."
    ;;
  n|N )
    echo "Operation aborted. PNG files were kept."
    exit 0
    ;;
  * )
    echo "Invalid input. No PNG files were deleted."
    exit 1
    ;;
esac
