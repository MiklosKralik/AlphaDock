#! /bin/bash
cmd=${1:-/bin/bash}
shift

#cd ~/AlphaDock
host=$(hostname|cut -f1 -d.)"-docker"
printf "\033]0;%s\007" docker

[ -z "$(pgrep docker)" ] && sudo service docker start
docker build -t alphadock ./docker_config
#container_id=$(docker run --rm -it -d -v $(pwd):/my -p 8080:80 -p 2222:22 -p 8888:8888 pandas-container)
container_id=$(docker ps | grep alphadock| head -1 | cut -f1 -d' ') 
if [ $container_id ]; then
    echo containers already running 
    docker ps
    echo attaching to $container_id
else
    container_id=$(docker run --rm -it -d -v $(pwd):/my -p 8888:8888 --hostname $host alphadock)
fi
echo -ne "\033]0;${USER}@${HOSTNAME}: ${PWD}\007"
docker exec -it $container_id $cmd $*
