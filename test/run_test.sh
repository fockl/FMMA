#!/bin/bash

blas_flag=0
if [[ $# > 0 ]]
then
  blas_flag=$1
fi

for name in ./build/*
do
  if [[ "$name" == *".out" ]]
  then
    printf "\e[34;1m ${name} \e[m\n"
    $name ${blas_flag}
    ret_code=$?
    if [ $ret_code != 0 ]
    then
      exit $ret_code
    fi
  fi
done
