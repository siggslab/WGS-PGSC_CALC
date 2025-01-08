#!/bin/bash


trace_file=$(ls -t results/runInfo/trace-*.tsv 2>/dev/null | head -n 1)

# Check if a trace file was found
if [ -z "$trace_file" ]; then
    echo "No trace-[timestamp].tsv files found in results/runInfo/"
    exit 1
fi

# Process the most recent trace file
#printf "$trace_file:\t"
# Find the header column index of 'workdir' and print the <workdir>/.command.log files
# Then, extract the SUs from the log files and sum them up
awk -v FS="\t" 'NR==1 { for (i=1; i<=NF; i++) if ($i=="workdir") wd=i } NR>1 { print $wd }' "$trace_file" | while read -r WORKDIR; do
    LOG="$WORKDIR/.command.log"
    if [ -f "$LOG" ]; then
        grep -m 1 "Service Units:" "$LOG" | sed -E -e 's/^.*:\s+([0-9\.]+)$/\1/'
    fi
done | awk -v FS="\t" '{ sum+=$0 } END { print sum }'
