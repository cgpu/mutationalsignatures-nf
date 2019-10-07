FROM nfcore/base:1.6
LABEL description="Docker image containing base requirements for mutational signatures analysis. Includes maftools,tidyverse, DT, ggplot2 and plotly."

COPY environment.yml /
RUN  conda env create -f /environment.yml && conda clean -a && \ 
     find /opt/conda/envs/mutationsignatures-nf/bin -type f -iname "*.R" -exec chmod +x {} \; && \
     find /opt/conda/envs/mutationsignatures-nf/bin -type f -iname "*.Rmd" -exec chmod +x {} \;

ENV PATH /opt/conda/envs/mutationsignatures-nf/bin:$PATH

