#!/bin/bash

mode=$(cat ./.mode)
sendfile=""

if test "$mode" = "go" ; then
  sendfile="main.go"
fi

if test "$mode" = "cpp"; then
  sendfile="main.cpp"
fi

dir=$(cat ./contest)
echo -e CONTEST: $dir
echo -n QUESTION?
read q
mkdir -p _result/_$dir/${q// /_}
mv -i $sendfile _result/_$dir/${q// /_}
git add -f _result/_$dir/${q// /_}/$sendfile
git commit -m "$dir $q"
cp -i _template/$sendfile ./$sendfile

