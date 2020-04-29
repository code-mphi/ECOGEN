#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# deploy only under travis request
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
    echo "Skipping deploy."
    exit 0
fi

#Clone the existing depot
echo "ECOGEN ----> Cloning repo ... "
git clone $GIT_DEPLOY_REPO temp
cd temp
git checkout $TARGET_BRANCH

#Cleaning old staff
rm -rf docs/sphinx_docs/*
rm -rf docs/doxygen_docs/*

#Building docs
echo "ECOGEN ----> Building Sphinx documentation ... "
cd ../docs/sphinx_docs
make html
cd ../..

echo "ECOGEN ----> Building Doxygen API documentation ... "
cd docs/doxygen_docs
doxygen doxygen.conf
cd ../..

#moving build in the good place
echo "ECOGEN ----> Moving build ... "
mkdir -p temp/docs/sphinx_docs
mv docs/sphinx_docs/build/html/* temp/docs/sphinx_docs
mkdir -p temp/docs/doxygen_docs
mv docs/doxygen_docs/build/html/* temp/docs/doxygen_docs

#committing changes
echo "ECOGEN ----> Commiting changes and push ... "
cd temp
git config user.name "Travis CI"
git config user.email "ecogen@code-mphi.fr"
git add .
git commit -m "Deploy to gh-pages"
git push $GIT_DEPLOY_REPO $TARGET_BRANCH:$TARGET_BRANCH
