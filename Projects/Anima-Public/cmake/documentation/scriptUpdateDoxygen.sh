#!/bin/bash

if [ ! -e ~/anima-build ]; then
	break
fi

cd ~/anima-build
make doc

nbFiles=`ls -l ~/anima-build/doc/html/* | wc -l`

if [ $nbFiles -le 2 ]; then
	break
fi

cd ~/dox-repo
rm -fr *
cp -r ~/anima-build/doc/html/* .

numChanges=`git diff | wc -l`

if [ $numChanges -eq 0 ]; then
	break
fi

echo git add --all
echo git commit -am "Cdash doxygen update"
echo git push origin gh-pages

