#!/bin/bash

for name in ./build/*
do
  if [[ "$name" == *".out" ]]
  then
    printf "\e[34;1m ${name} \e[m\n"
    $name
    ret_code=$?
    if [ $ret_code != 0 ]
    then
      exit $ret_code
    fi
  fi
done

echo "`cat comment_1.md`" > comment.md
echo "" >> comment.md
echo "` cat comment_2.md`" >> comment.md
