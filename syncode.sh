#! /bin/bash
# Copyright (C) 2016 jzzhao <jzzhao@zhlap>
# Distributed under terms of the MIT license.

target=u238:workspace/tcd
echo "Syncing $1 to $target ..."
case "$1" in
'cpp')
  rsync -azv CMakeLists.txt src $target/cpp
  ;;
'sh')
  rsync -azv --exclude 'syncode.sh' --delete-excluded *.sh $target/cpp
  ;;
'py')
  rsync -azv --exclude 'flycheck_*' --delete-excluded ../py/*.py $target/py
  ;;
*)
  echo "choose cpp, sh, or py"
  ;;
esac
