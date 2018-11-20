#!/usr/bin/env bash

read -r -d '' helpstring << EOM
Usage: docker run <merlot_container> [-h|--help] <mode> [mode options]
- scaffold: calculate a scaffold tree from manifold coordinates.
- elastic: calculate an elastic tree from a scaffold tree.
- embedded: embed an elastic tree back to gene space.
- inflate: inflate an elastic tree that was calculated from a reduced scaffold tree.
- -h|--help: print this message and exit.
EOM

key="$1"

case $key in
    scaffold)
    script="run_scaffold.R"
    ;;
    elastic)
    script="run_elastic.R"
    ;;
    embedded)
    script="run_embed.R"
    ;;
    inflate)
    script="inflate_elastic.R"
    ;;
    -h|--help)
    echo "$helpstring"
    ;;
    *)    # unknown option
    echo "Unknown mode. Please try again."
    echo "$helpstring"
    exit
    ;;
esac

Rscript /app/"$script" "${@:2}"
# # parse arguments again to figure out what to copy in and out
# POSITIONAL=()
# while [[ $# -gt 0 ]]
# do
#     key="$1"
#     echo "$key"
#     case $key in
#         -i|--input)
#         INFILE="$2"
#         shift # pass argument
#         shift # pass value
#         ;;
#         -o|--out)
#         OUTFILE="$2"
#         shift # pass argument
#         shift # pass value
#         ;;
#         -f|--full)
#         FULL="$2"
#         shift # pass argument
#         shift # pass value
#         ;;
#         *)    # unknown option
#         POSITIONAL+=("$1") # save it in an array for later
#         shift # pass argument
#         ;;
#     esac
# done

# set -- "${POSITIONAL[@]}"

# echo "$INFILE" "$OUTFILE" "$FULL"

# docker cp $INFILE tmp
# "${POSITIONAL[@]}"