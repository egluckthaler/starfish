#!/bin/bash

# inspired by https://nebulousresearch.org/other/bashcompletion

_starfish() {
	local cur
	_init_completion || return
	COMPREPLY=()
	cur="${COMP_WORDS[COMP_CWORD]}"
	opts="annotate consolidate sketch augment dereplicate-hood insert dereplicate flank extend summarize sim group cargo pair-viz locus-viz genome-viz format format-ncbi"
	COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
}
complete -F _starfish starfish

 