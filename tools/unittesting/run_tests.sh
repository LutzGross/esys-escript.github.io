#!/bin/bash

echo Pruning docker [ctrl-C to cancel].....
for i in {5..0} ;
do
	sleep 1
	echo -e "\e[1;31m $i \e[0m"
done
docker rmi `docker images -q` --force
docker rm `docker ps -q` --force
docker system prune --force
docker image prune --force
docker container prune --force
docker volume prune --force

echo -e "\e[1;31m  Starting testing..... \e[0m"
echo -e "\e[1;31m  Testing debian PASO ..... \e[0m"
docker build -f debian_py3 .          | tee out.debian_py3_test
echo -e "\e[1;31m  Testing debian Trilinos ..... \e[0m"
docker build -f debian_trilinos_py3 . | tee out.debian_trilinos_py3_test
# echo -e "\e[1;31m  Testing opensuse PASO ..... \e[0m"
# docker build -f opensuse_py3 .          | tee out.opensuse_py3_test
# echo -e "\e[1;31m  Testing opensuse Trilinos ..... \e[0m"
# docker build -f opensuse_trilinos_py3 . | tee out.opensuse_trilinos_py3_test
echo -e "\e[1;31m  Testing opensuse MPI ..... \e[0m"
docker build -f opensuse_trilinos_MPI_py3 . | tee out.opensuse_trilinos_MPI_py3_test

