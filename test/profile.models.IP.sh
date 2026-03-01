#!/bin/zsh

start="${1:-3}"
stop="${2:-4}"
dim="${3:-4}"
fields=(2)

echo "________________________________"
echo "| PROFILE INVASION PERCOLATION | ➭➭➭ results in profiles/InvasionPercolation"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

source .vars

printf "%-5s %-${W}s %-${W}s %-${W}s %-${W}s %s" "" $SCALE $DIM $FIELD $PHAT $CELLS
echo

for ((L=$start; L<$stop; L++)); do
	python profile.models.IP.py $dim $L 2 &
	wait $!
	
	case $? in
		1)
			echo -e "$UP$FAIL";;
		0)
			echo -e "$UP$PASS";;
		*)
			echo -e "$UP$WARN";;
	esac
done

width=$(tput cols)
indent=2

echo
echo