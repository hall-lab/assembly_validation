FROM rocker/tidyverse@sha256:0d900477d9ca90ff297db9f444a5f8ea641f70463524d5804a46d89dc141093a
LABEL maintainer "Allison Regier <aregier@wustl.edu>"

RUN mkdir -p /opt/hall-lab
ADD make_plots.R /opt/hall-lab

CMD ["/bin/bash"]
