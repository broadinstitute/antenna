#!/bin/bash

tag=0.0.10

if [ -z $tag ]; then
    echo -e "\nYou must provide a tag"
    echo -e "\nUsage: bash build_docker.sh TAG\n"
    exit 1
fi

docker build -t quay.io/nbarkas_1/antenna:$tag .
