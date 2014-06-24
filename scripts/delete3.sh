#!/bin/bash

for f in AC*
do
    sed -i -e 1,3d $f
done
