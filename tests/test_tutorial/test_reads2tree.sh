#!/usr/bin/env bash
parallel=true
install=true
while [ "$#" -gt 0 ]; do
    case "$1" in
	"--no-parallel" | "-P")
	    parallel=false
	    ;;
	"--no-install" | "-I")
	    install=false
	    ;;
	"--assemblies" | "-a")
	    shift
	    assemblies_dir="$(realpath "$1")"
	    ;;
	*)
	    >&2 echo "Unrecognized argument $1."
	    exit 1
	    ;;
    esac
    shift
done
if ! which sudo; then
    function sudo {
	"$@"
    }
fi
cd
set -e
case "$(uname)" in
    Linux)
	if [ "$parallel" = true ] && [ "$install" = true ]; then
	    sudo apt update
	    sudo NEEDRESTART_MODE=a apt -y install parallel
	fi
	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
	;;
    Darwin)
	if [ "$parallel" = true ] && [ "$install" = true ]; then
	    brew install parallel
	fi
	case "$(uname -m)" in
	    x86_64)
		curl -L -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
		;;
	    arm64)
		curl -L -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac-arm64.tar.gz
		;;
	    *)
		>&2 echo "Unrecognized machine type $(uname -m)."
		exit 1
		;;
	esac
	;;
    *)
	>&2 echo "Unrecognized system $(uname)."
	exit 1
	;;
esac
tar xzvf sratoolkit.current-*.tar.gz
export PATH="$PATH:$(realpath sratoolkit*/bin)"
eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
conda activate rna-clique
conda install -y lxml requests
git clone https://github.com/actapia/download_sra
export PATH="$PATH:$PWD/download_sra"
case "$(uname)" in
    Linux)
	wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
	;;
    Darwin)
	curl -L -O https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Darwin-$(uname -m).tar.gz
	;;
    *)
	>&2 echo "Unrecognized system $(uname)."
	exit 1
	;;
esac
tar xzvf SPAdes-*.tar.gz
PATH="$PATH:$(realpath SPAdes-*/bin)"
export PATH
cd rna_clique
export RNA_CLIQUE="$PWD"
cd ..
mkdir tutorial
cd tutorial
export TUTORIAL_DIR="$PWD"
cd "$TUTORIAL_DIR"
cut -d, -f1 "$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv" | \
    download_sra.sh -j 0 -r
compgen -G "SRR*.fastq"
if [[ -v assemblies_dir ]]; then
    ln -s "$assemblies_dir" "out"
else
    if [ "$parallel" = true ]; then
	parallel --jobs 1 spades.py --rna -o out/{/.} -s {} -t 3 -m 13 ::: *.fastq
    else
	for f in *.fastq; do
	    b="$(basename "$f")"; 
	    fn="${b%%.*}";
	    spades.py --rna -o "out/$fn" -s "$f" -t 3 -m 120;
	done
    fi
fi
for f in out/*; do
    [ -f "$f/transcripts.fasta" ]
done
cd "$RNA_CLIQUE"
out_dir="$TUTORIAL_DIR/rna_clique_out/"
bash typical_filtering_step.sh -o "$out_dir" \
     -n 50000 \
     "$TUTORIAL_DIR"/out/*
[ -d "$out_dir" ]
[ -f "$out_dir/graph.pkl" ]
PYTHONPATH='.' python docs/tutorials/reads2tree/make_tree.py
[ -f "$out_dir/nj_tree.svg" ]
[ -f "$out_dir/nj_tree.tree" ]
PYTHONPATH="." python docs/tutorials/reads2tree/make_pcoa.py
[ -f "$out_dir/pcoa_2d.svg" ]
[ -f "$out_dir/pcoa_3d.svg" ]
PYTHONPATH="." python docs/tutorials/reads2tree/make_heatmap.py
[ -f "$out_dir/distance_heatmap.svg" ]

