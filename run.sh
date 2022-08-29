#!/bin/bash

docker run --rm --name gravity2 -p 8000:80 \
    -e PATH=$PATH \
    -e HHLIB=$HHLIB \
    gravity2/