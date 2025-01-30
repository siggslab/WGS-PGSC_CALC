#!/bin/bash



trace_file=$(ls -t $1 2>/dev/null | head -n 1)

# Check if a trace file was found
if [ -z "$trace_file" ]; then
    echo "No trace-[timestamp].tsv files found in results/runInfo/"
    exit 1
fi

# Process the most recent trace file
#printf "$trace_file:\t"
# Find the header column index of 'workdir' and print the <workdir>/.command.log files
# Then, extract the SUs from the log files and sum them up
# awk -v FS="\t" 'NR==1 { for (i=1; i<=NF; i++) if ($i=="workdir") wd=i } NR>1 { print $wd }' "$trace_file" | while read -r WORKDIR; do
#     LOG="$WORKDIR/.command.log"
#     if [ -f "$LOG" ]; then
#         grep -m 1 "Service Units:" "$LOG" | sed -E -e 's/^.*:\s+([0-9\.]+)$/\1/'
#     fi
# done | awk -v FS="\t" '{ sum+=$0 } END { print sum }'
# Initialize sums for cached and non-cached SUs


# Process the trace file using awk
awk -v FS="\t" '
    BEGIN {
        cached_sum = 0
        non_cached_sum = 0
    }
    NR == 1 {
        # Identify column indices for status and workdir
        for (i = 1; i <= NF; i++) {
            if ($i == "status") status_idx = i
            if ($i == "workdir") workdir_idx = i
        }
    }
    NR > 1 {
        # Construct the path to .command.log
        log_file = $workdir_idx "/.command.log"

        # Check if the log file exists and extract SU value
        if ((cmd = "grep -m 1 \"Service Units:\" " log_file " 2>/dev/null") | getline line) {
            close(cmd)
            if (match(line, /[0-9.]+/)) {
                SU = substr(line, RSTART, RLENGTH)

                # Add SU to the appropriate sum based on status
                if ($status_idx == "CACHED") {
                    cached_sum += SU
                } else {
                    non_cached_sum += SU
                }
            }
        }
    }
    END {
        # Print the results
        total_sum = cached_sum + non_cached_sum
        print total_sum " \033[90m(" non_cached_sum " new, " cached_sum " cached)\033[0m"
    }
' "$trace_file"



# Print the results
# echo "Cached SU usage: $cached_sum"
# echo "Non-cached SU usage: $non_cached_sum"
# echo "$(echo "$cached_sum + $non_cached_sum" | bc) ($cached_sum cached)"