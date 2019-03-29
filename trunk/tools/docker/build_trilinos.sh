#!/bin/bash

if [ "$EUID" -ne 0 ]
	then echo Please run as sudo
	exit
fi

docker build --tag=escript5.3 -f docker_trilinos .
