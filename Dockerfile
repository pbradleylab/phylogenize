FROM condaforge/mambaforge
LABEL org.opencontainers.image.authors="kananen.13@osu.edu"

RUN mamba install -y -c bioconda \
	bioconductor-qvalue \
	bioconductor-ggtree \
	bioconductor-biomformat \
	vsearch
RUN mamba install -y -c conda-forge \
	r-devtools \
	r-ragg \
	r-phylolm \
	r-phangorn

#RUN R -e "devtools::install_bitbucket('pbradz/phylogenize/package/phylogenize')"
RUN R -e "devtools::install_github('pbradleylab/phylogenize')"

CMD ["/bin/bash"]
