#!/bin/bash
# rename_prefix.sh
# Rename all files starting with prefix OLD to prefix NEW in a directory

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 OLD_PREFIX NEW_PREFIX [DIRECTORY]"
  exit 1
fi

OLD="$1"
NEW="$2"
DIR="${3:-.}"   # default: current directory

# Safety check
if [ ! -d "$DIR" ]; then
  echo "Error: $DIR is not a directory"
  exit 1
fi

shopt -s nullglob
for f in "$DIR"/"$OLD"*; do
  [ -f "$f" ] || continue
  base="$(basename "$f")"
  new="$NEW${base#$OLD}"
  mv -- "$f" "$DIR/$new"
  echo "Renamed: $base → $new"
done
