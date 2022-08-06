#! awk -f
# converts DNA strings with CIGAR strings into aligned strings
#
# expected input, lines of:
# [sequence]	[cigar string]
# allowing for '-' as cigar of reference sequence

# transfor '- cigar into X repeat
function expand_ref_cigar() {
	for (i in cigs) {
		if (cigs[i]=="-") {
			cigs[i] = length(seqs[i]) "X"
		}
	}
}

# parse vector of cigar strings into arrays com/cts
# cigs: cigar strings
# com: commands vector, with M, X, D, I
# cts: counts vector, with integers
# lengths : size of com[i *]/cts[i *]
function break_cigar() {
	for (i in cigs) {
		split(cigs[i], comi, /[0-9]+/, ctsi)
		delete comi[1]
		# storing elements in com[i j] and cts[i j]
		lengths[i] = length(comi)
		for (j=1; j<=length(comi); j++) {
			com[i j] = comi[j+1]
			cts[i j] = ctsi[j]
		}
	}
}

# print current state of cigar strings, in arrays com/cts
function print_cigar() {
	for(i=1; i<=N; i++) {
		for(j=1; j<=lengths[i]; j++) {
			printf "%d%s,", cts[i j], com[i j]
		}
		printf "\n"
	}
}

# print vector, one element at a line
function print_array(ar,    i) {
	for(i in ar) {
		printf "%s\n", ar[i]
	}
}

# rep(char, n) returns string of char repeated n times
function rep(char, n) {
	ret=""
	for (i=1; i<=n; i++) ret = ret char
	return(ret)
}

# find current (at icig) I commands, returning in v[] vector of their
# indexes, but in order of their cts size (in case more than one present)
function current_Is(v, icig,     k, savei) {
	v[0] = 0 # number of I commands found
	for (k=1; k<=N; k++) {
		if (com[k icig[k]]=="I") {
			# inserting in increasing order
			if (v[0]==0) v[1] = k
			else {
				if (cts[k icig[k]] < cts[v[v[0]] icig[v[v[0]]]]) {
					savei = v[v[0]]
					v[v[0]] = k
					v[v[0]+1] = savei
				} else {
					v[v[0]+1] = k
				}
			}
			v[0]++
		}
	}
}

## strategy: to parse cigar strings one character at a time
function naive_parser(  ended, i, done, v, toadd, out, iseq, icig) {
	# initializing vectors
	for (i=1; i<=N; i++) {
		out[i]=""
		iseq[i]=1
		icig[i]=1
	}
	ended = 0
	#while (ended<length(seqs[1])) {
	while (!ended) {
		# first, deals with I commands
		v[0] = 0
		current_Is(v, icig)
		if (v[0]>0) {
			done = 0 # how many were inserted already
			for (i=1; i<=v[0]; i++) {
				# process command I at cigs[v[i] icig[i]]
				toadd = cts[v[i] icig[v[i]]] - done
				done = toadd + done
				for (j=1; j<=N; j++) {
					if (j!=v[i]) {
						out[j] = sprintf("%s%s", out[j], rep(" ", toadd))
					} else {
						out[j] = sprintf("%s%s", out[j], substr(seqs[j], iseq[j], cts[j icig[j]]))
						iseq[j]+=cts[j icig[j]]
					}
				}
				cts[v[i] icig[v[i]]] = 0
				icig[v[i]]++
			}
		}
		# other commands
		for (i=1; i<=N; i++) {
			if (com[i icig[i]] == "X") {
				out[i] = sprintf("%s%s", out[i], substr(seqs[i], iseq[i], 1))
				iseq[i]++
			} else if (com[i icig[i]] == "M") {
				out[i] = sprintf("%s-", out[i])
				iseq[i]++
			} else if (com[i icig[i]] == "D") {
				out[i] = sprintf("%s ", out[i])
			}
			else { print "--error--" }
			cts[i icig[i]]--
			if (cts[i icig[i]]==0) icig[i]++
			# check processed all commands
			if (icig[i]>lengths[i]) ended = 1
		}
	}
	# result stored back to seqs array
	for (i=1; i<=N; i++) {
		seqs[i] = out[i]
	}
}

#### main()
BEGIN{N=0} # N: number of cigar strings
{
	N++
	seqs[N] = $1
	cigs[N] = $2
}
END{
	expand_ref_cigar()
	break_cigar()
	# printf "initial state:\n"
	# print_cigar()
	# print_array(seqs)
	naive_parser()
	
	# printf "final state:\n"
	# print_cigar()
	# printf "final sequences:\n"
	print_array(seqs)
}
