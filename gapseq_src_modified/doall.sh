#!/bin/bash

bitcutoff=200
min_bs_core=100
identcutoff=0   # cutoff blast: min identity
covcutoff=75 # cutoff blast: min coverage
output_dir=.
user_temp=false
aliTool="blast"
aliTool_user_set=false
medium="auto"
taxonomy="auto"
user_temp_opt=""
verbose=0
force_offline=false
cache_dir=${GAPSEQ_CACHE_DIR:-$HOME/.gapseq_cache}
cluster_mode="off"
target_cache_dir=${GAPSEQ_TARGET_CACHE_DIR:-$HOME/.gapseq_target_cache}
target_cache_mode="off"
ecoli_fast_db=""

userdir=false

OS=$(uname -s)
if [ "$OS" = "Darwin" -o "$OS" = "FreeBSD" ]; then
	n_threads=$(sysctl hw.ncpu|cut -f2 -d' ')
else
	n_threads=`grep -c ^processor /proc/cpuinfo`
fi

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
seqdb=$dir/../dat/seq

usage()
{
    echo "gapseq - perfoming all reconstruction steps (doall)"
    echo ""
    echo "Usage"
    echo "$0 file.fasta"
    echo "  -b bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -l Reactions with an associated blast-hit with a bitscore below this value will be considered just as reactions that have no blast hit."
    echo "  -t Taxonomic range for reference sequences to be used. (Bacteria, Archaea, auto; default: $taxonomy)."
    echo "  -m tab- or comma-separated table for media components. Requires three named columns: 1 - \"compounds\" (for metab. IDs), 2 - \"name\" (metab. name), 3 - \"maxFlux\" (maximum inflow flux). If \"auto\" (the default), a growth medium will be predicted."
    echo "  -f Path to directory, where output files will be saved (default: current directory)"
    echo "  -v Verbose level, 0 for nothing, 1 for full (default $verbose)"
    echo "  -D path to directory, where reference sequence database will be saved (default: $seqdb)"
    echo "  -T Set user-defined temporary folder (default: $user_temp)"
    echo "  -O Force offline mode (default: $force_offline)"
    echo "  -K Number of threads. If option is not provided, number of available CPUs will be automatically determined."
    echo "  -A Tool to be used for sequence alignments (blast, mmseqs2, diamond; default: $aliTool)"
    echo "  -C Cache directory for alignment results (default: $cache_dir)"
    echo "  -S Cluster identical query proteins before alignment (off|identical; default: $cluster_mode)"
    echo "  -G Enable cross-genome target cache (off|on) and optionally set cache dir (default: $target_cache_mode)"
    echo "  -E E. coli pan-reference FASTA for fast+rescue alignment mode in gapseq_find"
    echo ""
exit 1
}

OPTIND=1         # Reset in case getopts has been used previously in the shell.
while getopts "h?b:i:c:l:t:m:f:v:OD:T:K:A:C:S:G:E:" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    b)
        bitcutoff=$OPTARG
        ;;
    i)
        identcutoff=$OPTARG
        ;;
    c)
        covcutoff=$OPTARG
        ;;
    l)
        min_bs_core=$OPTARG
        ;;
    t)
        taxonomy=$OPTARG
        ;;
    m)
        medium=$OPTARG
        ;;
    f)
        output_dir=$OPTARG
        ;;
    v)
        verbose=$OPTARG
        ;;
    D)
        seqdb=$OPTARG
        userdir=true
        ;;
    T)
        user_temp=true
        user_temp_folder=$OPTARG
        user_temp_opt="-T ${user_temp_folder}"
        ;;
    O)
        force_offline=true
        ;;
    K)
        n_threads=$OPTARG
        ;;
    A)
        aliTool=$(echo "$OPTARG" | tr '[:upper:]' '[:lower:]')
        aliTool_user_set=true
        if [[ "$aliTool" != "blast" && "$aliTool" != "diamond" && "$aliTool" != "mmseqs2" && "$aliTool" != "tblastn" ]]; then
            echo "Error: Invalid value for -A. Expected 'blast', 'diamond', 'mmseqs2', or 'tblastn', got '$OPTARG'."
            exit 1
        fi
        ;;
    C)
        cache_dir="$OPTARG"
        ;;
    S)
        cluster_mode=$(echo "$OPTARG" | tr '[:upper:]' '[:lower:]')
        if [[ "$cluster_mode" != "off" && "$cluster_mode" != "identical" ]]; then
            echo "Error: Invalid value for -S. Expected 'off' or 'identical', got '$OPTARG'."
            exit 1
        fi
        ;;
    G)
        if [[ "$OPTARG" == "on" || "$OPTARG" == "off" ]]; then
            target_cache_mode="$OPTARG"
        else
            target_cache_mode="on"
            target_cache_dir="$OPTARG"
        fi
        ;;
    E)
        ecoli_fast_db=$(readlink -f "$OPTARG")
        ;;
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift
# after parsing arguments, only fasta file shoud be there
[ "$#" -ne 1 ] && { usage; }

# Auto-select fastest available aligner when user hasn't explicitly set -A
if [ "$aliTool_user_set" = false ]; then
    if command -v diamond >/dev/null 2>&1; then
        aliTool="diamond"
        echo "Auto-selected diamond for faster alignment (use -A blast to override)."
    fi
fi

# Ensure selected aligner is available before starting heavy computation.
if [ "$aliTool" = "blast" ] && ! command -v blastp >/dev/null 2>&1; then
    echo "Error: blastp not found in PATH. Please install BLAST+ or set -A to another available aligner."
    exit 1
fi
if [ "$aliTool" = "diamond" ] && ! command -v diamond >/dev/null 2>&1; then
    echo "Error: diamond not found in PATH. Please install DIAMOND or set -A blast."
    exit 1
fi


# set output directory
case $output_dir in
    /*)
        # echo "absolute path"
        output_dir=$output_dir
        ;;
    ~*)
        # echo "relative to $HOME directory"
        output_dir="${output_dir/#\~/$HOME}"
        ;;
    *)
        # echo "relative path to current directory"
        output_dir=$curdir/$output_dir
        ;;
esac

file=$(readlink -f $1)
base=$(basename "$file")
id="${base%.*}"
[[ ! -s "$file" ]]  && usage
[[ $file == *.gz ]] && id="${id%.*}"

# Avoid CPU oversubscription when running both stages in parallel.
find_threads=$n_threads
trans_threads=$n_threads
if [ "$n_threads" -ge 2 ]; then
    find_threads=$(( (n_threads + 1) / 2 ))
    trans_threads=$(( n_threads / 2 ))
    [ "$trans_threads" -lt 1 ] && trans_threads=1
fi

# Run gapseq_find and transporter in parallel for performance
$dir/gapseq_find.sh -v $verbose -b $bitcutoff -p all -t $taxonomy -K $find_threads -A $aliTool -C "$cache_dir" -S "$cluster_mode" -G "$target_cache_mode" -f $output_dir $user_temp_opt $([ "$force_offline" = true ] && echo -O) $([ "$userdir" = true ] && echo -D $seqdb) $([ -n "$ecoli_fast_db" ] && echo "-E $ecoli_fast_db") "$file" &
pid_find=$!

$dir/transporter.sh -v $verbose -b $bitcutoff -K $trans_threads -A $aliTool -C "$cache_dir" -S "$cluster_mode" -G "$target_cache_mode" -f $output_dir $user_temp_opt "$file" &
pid_trans=$!

# Wait for both to complete
wait $pid_find
exit_find=$?
wait $pid_trans
exit_trans=$?

# Check if either failed
if [ $exit_find -ne 0 ]; then
    echo "Error: gapseq_find.sh failed with exit code $exit_find"
    exit $exit_find
fi
if [ $exit_trans -ne 0 ]; then
    echo "Error: transporter.sh failed with exit code $exit_trans"
    exit $exit_trans
fi

Rscript $dir/generate_GSdraft.R -r "$output_dir/$id-all-Reactions.tbl" -t "$output_dir/$id-Transporter.tbl" -c "$file" -u $bitcutoff -l $min_bs_core -p "$output_dir/$id-all-Pathways.tbl" -b $taxonomy -f $output_dir
if [ $medium == "auto" ]; then
    Rscript $dir/predict_medium.R -m "$output_dir/${id}-draft.RDS" -p "$output_dir/$id-all-Pathways.tbl" -f $output_dir
    medium="$output_dir/${id}-medium.csv"
fi
Rscript $dir/gf.suite.R -m "$output_dir/${id}-draft.RDS" -n "$medium" -c "$output_dir/${id}-rxnWeights.RDS" -g "$output_dir/${id}-rxnXgenes.RDS" -b $min_bs_core -f $output_dir


