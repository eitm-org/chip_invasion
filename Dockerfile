FROM rocker/shiny:latest

LABEL maintainer="Elizabeth Elton <eelton@eitm.org>"

# system libraries of general use
RUN apt-get update && apt-get install --no-install-recommends -y \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.1 \
    libxml2 \
    && rm -rf /var/lib/apt/lists/*

# system library dependency for the example app
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

# basic shiny functionality
RUN R -q -e "install.packages(c('shiny', 'rmarkdown', 'bslib'))"

# install dependencies of the example app
RUN R -q -e "install.packages(c('tidyverse', 'here', 'janitor', 'readxl', 'hablar', 'ggbeeswarm', 'plotly', 'FactoMineR', 'factoextra','xml2'))"

# copy the app to the image
RUN mkdir /root/chip_invasion
COPY chip_invasion /root/chip_invasion

COPY Rprofile.site /usr/local/lib/R/etc/

EXPOSE 3838

# uncomment this line to run on bare shiny instead of shinyproxy
# CMD ["R", "-q", "-e", "shiny::runApp('/root/chip_invasion')"]
