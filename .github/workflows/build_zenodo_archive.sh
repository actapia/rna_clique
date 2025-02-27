#!/usr/bin/env bash
set -x
set -e
set -o pipefail
declare -a urls
git clone https://github.com/actapia/download_sra/
python git_archive_all_extra.py "$1.zip"
while IFS= read -r -u 5 line; do
    urls+=("$(python download_sra/get_run_link.py "$line" | tail -n1)")
done 5< <(tail -n+2 "$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv" | cut -d, -f1)
python add_remote_files_to_archive.py "$1.zip" --remote "${urls[@]}" --remote-dir "$1/test_data/sra"
python add_remote_files_to_archive.py "$1.zip" --remote "http://rna-clique-data.s3-website.us-east-2.amazonaws.com/transcripts.zip" --remote-dir "$1/test_data/assemblies" --merge-zips --compress-level 9
mv "$1.zip" "$RNA_CLIQUE"
