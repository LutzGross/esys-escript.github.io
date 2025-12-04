#!/bin/bash
# Script to build and publish full documentation to GitHub Pages

set -e

echo "Building full documentation..."
scons -j4 docs

echo "Checking out gh-pages branch..."
git fetch origin
if git rev-parse --verify gh-pages >/dev/null 2>&1; then
    git branch -D gh-pages || true
fi
git checkout --orphan gh-pages

echo "Cleaning working directory..."
git rm -rf .
git clean -fdx

echo "Copying documentation..."
cp -r release/doc/* .

echo "Adding .nojekyll file..."
touch .nojekyll

echo "Committing documentation..."
git add .
git commit -m "Deploy documentation

Built from commit: $(git rev-parse master)
Date: $(date -u '+%Y-%m-%d %H:%M:%S UTC')

🤖 Generated with [Claude Code](https://claude.com/claude-code)
"

echo "Pushing to gh-pages branch..."
git push origin gh-pages --force

echo "Returning to master branch..."
git checkout master

echo ""
echo "Documentation published successfully!"
echo "Visit: https://lutzgross.github.io/esys-escript.github.io/"
