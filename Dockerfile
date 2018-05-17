FROM rocker/shiny:latest

RUN apt-get update -y && \
    apt-get install libssl-dev libxml2-dev libnetcdf-dev -y && \
    echo "install.packages('igraph')\ninstall.packages('reshape')" | r && \
    echo "library(utils)\nsource('http://bioconductor.org/biocLite.R')\nbiocLite('SIMLR')\nbiocLite('scater')\nbiocLite('monocle')" | r

WORKDIR /srv/shiny-server/SingleCell

ADD R /srv/shiny-server/SingleCell

RUN chown -R shiny:shiny . && mkdir data

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]