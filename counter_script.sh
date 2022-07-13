#!/bin/bash

START_DIR=$(pwd)

for dir in */
do
cd $dir
pwd
find . | wc -l
find . -type f | wc -l
cd $START_DIR
done

