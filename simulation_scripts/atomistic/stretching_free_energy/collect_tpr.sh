#!/bin/bash

start_int=0
end_int=1540
step_int=2

# Loop over coord_* directories
current=$start_int
while (( $(echo "$current <= $end_int" | bc -l) )); do
    # Format directory name
    #dir=$(printf "coord_%.4f" "$current")
    dir=coord_$current
    # Check if directory exists
    if [ -d "$dir" ]; then
        # Find pro*.tpr and pro*_pullf.xvg files
        for pullf in "$dir"/OUT*_pullf.xvg; do
            # Check if pro*_pullf.xvg is non-empty
            if [ -s "$pullf" ]; then
                # Find the corresponding .tpr file
                tpr="${pullf%_pullf.xvg}.tpr"

                # Check if the .tpr file exists
                if [ -f "$tpr" ]; then
                    # Append the names to respective files
                    echo "$tpr" >> tpr.dat
                    echo "$pullf" >> pullf.dat
                fi
            fi
        done
    fi

    # Increment the current value
    current=$(echo "$current + $step_int" | bc)
done



echo "Done. Results saved to tpr.dat and pullf.dat."

