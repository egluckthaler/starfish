#!/bin/bash

# inspired by https://nebulousresearch.org/other/bashcompletion
# and https://github.com/soedinglab/MMseqs2/blob/master/util/bash-completion.sh

_starfish() {
	local cur
	_init_completion || return
	COMPREPLY=()
	cur="${COMP_WORDS[COMP_CWORD]}"

	if [[ ${COMP_CWORD} -eq 1 ]] ; then
		opts="annotate consolidate coverage sketch augment dereplicate-hood insert dereplicate flank extend summarize sim group cargo pair-viz locus-viz genome-viz format format-ncbi"
		COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
		return 0
	fi
	
	if [[ ${COMP_CWORD} -eq 2 ]] ; then
		opts="-h"
		COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
		return 0
	fi

	if [[ ${COMP_CWORD} -gt 2 ]] ; then
		compopt -o default
		COMPREPLY=( )
		return 0
	fi
	
}
complete -F _starfish starfish

 