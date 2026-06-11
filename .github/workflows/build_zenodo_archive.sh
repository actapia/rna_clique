#!/usr/bin/env bash
set -x
set -e
set -o pipefail
declare -a urls
git clone https://github.com/actapia/download_sra/
python git_archive_all_extra.py "$1.zip"
[ -f "$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv" ]
while IFS= read -r -u 5 line; do
    urls+=("$(python download_sra/get_run_link.py "$line" | tail -n1)")
done 5< <(tail -n+2 "$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv" | cut -d, -f1)
python add_remote_files_to_archive.py "$1.zip" --remote "${urls[@]}" --remote-dir "$1/test_data/sra"
python add_remote_files_to_archive.py "$1.zip" --remote "http://rna-clique-data.s3-website.us-east-2.amazonaws.com/transcripts.tar.gz" --remote-dir "$1/test_data/assemblies" --merge-zips --compress-level 9
python add_remote_files_to_archive.py "$1.zip" --remote "http://rna-clique-data.s3-website.us-east-2.amazonaws.com/trinity_assemblies.tar.xz" --remote-dir "$1/test_data/trinity_assemblies" --merge-zips --compress-level 9
python add_remote_files_to_archive.py "$1.zip" --remote "https://cs.uky.edu/~acta225/rna_clique/ec_genomes.tar.xz" --remote-dir "$1/test_data/ec_genomes" --merge-zips --compress-level 9
mv "$1.zip" "$RNA_CLIQUE"
