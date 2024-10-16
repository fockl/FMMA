#!/bin/bash

for name in ./build/*
do
  echo $name
  if [[ "$name" == *".out" ]]
  then
    echo $name
    printf "\e[34;1m ${name} \e[m\n"
    $name
    ret_code=$?
    if [ $ret_code != 0 ]
    then
      exit $ret_code
    fi
  fi
done
