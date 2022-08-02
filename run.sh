#!/bin/bash

docker run --rm --name gravity2 -p 8000:80\
    -e ENV_VAR=$MAPPING\
    GRAViTy2/${PWD##*/}