#!/bin/bash

# script courtesy of T. Hackl. Downloaded from https://github.com/thackl/seq-scripts/blob/master/bin/seq-gc 29/09/22

## subs
usage(){
cat <<EOF
Usage:
  seq-gc FASTA > FASTA.gc.tsv
Calculate per sequence / window GC content.
  -w  Window size. Unless set, GC is calculated per sequence.
  -s  Shift between window starts. Default -w.
  -m  Calculate median of windows. Default -w=1000.
  -b  Output bed instead of tsv format for -w mode.
  -N  ignore Ns (and other non ATGC) in per sequence mode
      gives AT+GC = 100% even on scaffolds, 0 for N only.
EOF
exit 0;
}

check_bin(){
    hash $1 || { echo "$1 required in PATH" >&2; exit 1;}
}

## prep
[[ $# -eq 0 ]] && usage;

# Execute getopt
ARGS=`getopt --name "seq-gc" \
    --options "w:s:bNmhV" \
    -- "$@"`

#Bad arguments
[ $? -ne 0 ] && exit 1;

# A little magic
eval set -- "$ARGS"

# Now go through all the options
while true; do
    case "$1" in
        -w)
            [ ! -n "$2" ] && (echo "$1: value required" 1>&2 && exit 1);
            WIN=$2;
            shift 2;;
        -s)
            [ ! -n "$2" ] && (echo "$1: value required" 1>&2 && exit 1);
            SHIFT=$2;
            shift 2;;
        -m)
            MEDIAN=1;
            shift;;
        -b)
            BED=1;
            shift;;
        -N)
            NIGNORE=1;
            shift;;
        -h)
	    usage && exit 0;;

        -V)
            grep -m1 'Version' "$0" | sed 's/.*Version\s*//';
            exit 0;;
        --)
            shift
            break;;
        *)
            echo "$1: Unknown option" 1>&2 && exit 1;;
    esac
done


## prep
[[ $# -eq 0 ]] && usage;

check_bin samtools;
check_bin bedtools;

FA=$1;
PRE=`basename $FA`
PRE=.${PRE%.*};
[[ -n $MEDIAN ]] && [[ -z $WIN ]] && WIN=1000;
[[ -z $SHIFT ]] && SHIFT=$WIN

## main
samtools faidx $FA;

if [[ -z $WIN ]]; then
    cut -f 1,2 $FA.fai | sed 's/\t/\t0\t/' > $PRE-gc.bed
else
    cut -f 1,2 $FA.fai > $PRE.len
    bedtools makewindows \
        -g $PRE.len \
        -w $WIN \
        -s $SHIFT \
        > $PRE-gc.bed
fi;

if [[ -z $MEDIAN ]]; then  # no median
    if [[ -z $NIGNORE ]]; then  # dont ignore Ns
        if [[ -z $BED ]]; then  # output bed
            bedtools nuc \
                     -fi $FA \
                     -bed $PRE-gc.bed |
                cut -f1,5 |
                tail -n +2
        else
            bedtools nuc \
                     -fi $FA \
                     -bed $PRE-gc.bed |
        gawk -v w=${WIN} 'BEGIN{FS="\t"; OFS="\t"}
{
if (FNR>1) {print $1,$2,$3,"GC.w"w,$5}
}'
        fi;
    else  # ignore Ns
        if [[ -z $BED ]]; then  # output bed
        bedtools nuc \
            -fi $FA \
            -bed $PRE-gc.bed |
            tail -n +2 |
            perl -ane '
              my $gcc = $F[6] + $F[7];
              my $tot = $gcc + $F[5] + $F[8];
              printf "$F[0]\t%0.6f\n", $tot ? $gcc/$tot : 0;
            '
        else
            bedtools nuc \
                     -fi $FA \
                     -bed $PRE-gc.bed |
                tail -n +2 |
                perl -ane '
                  my $gcc = $F[6] + $F[7];
                  my $tot = $gcc + $F[5] + $F[8];
                  printf "%s\t%s\t%s\tGC.w'$WIN'\t%0.6f\n", @F[0..2], $tot ? $gcc/$tot : 0;
                '
        fi;
    fi;
else
    bedtools nuc \
        -fi $FA \
        -bed $PRE-gc.bed |
    cut -f1,5 |
    tail -n +2 |
    perl -ane '
      if($F[0] eq $id){
          push @w, $F[1]
      }else{
          if(@w){
              @w=sort{$a<=>$b}@w;
              print $id,"\t",$w[@w/2],"\n"
          }
          $id=$F[0];
          @w=($F[1]);
      };
      END{
          if(@w){
              @w=sort{$a<=>$b}@w;
              print $id,"\t",$w[@w/2],"\n"
          }
      }'
fi;

rm -f $PRE-gc.bed $PRE.len;