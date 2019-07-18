#!/bin/bash

if [ "$EUID" -ne 0 ]
	then echo Please run as sudo
	exit
fi

docker build --tag=esysescript/esys-escript -f docker_debian .
