#!/bin/bash

mkdir -p data/analysistree
rsync -av /u/klochkov/temp/temp/atree/1.analysistree.root data/analysistree/

ls -d data/analysistree/ > data/filelist.txt