#!/bin/bash

cd ./calculator
rm *
cp /home/oleksii/cbmdir/flow_calculator/macro/* ./

cd ../discriminator
rm *
cp /home/oleksii/cbmdir/qn_discriminator/macro/* ./

cd ../drawing_tools
rm *
cp /home/oleksii/cbmdir/flow_drawing_tools/macro/* ./

cd ..

git add calculator discriminator drawing_tools

if [[ ${1} = "" ]]
then
MESSAGE="(ir)regular update"
else
MESSAGE=${1}
fi

git commit -m "$MESSAGE"
git push
