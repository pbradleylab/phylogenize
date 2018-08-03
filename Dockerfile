FROM r-base
WORKDIR /usr/local/src
RUN git clone https://pbradz@bitbucket.org/pbradz/phylogenize.git
RUN wget https://github.com/knights-lab/BURST/releases/download/v0.99.7f/burst-0.99.7f-linux-64.tar.gz
RUN tar xvfz burst-0.99.7f*tar.gz
COPY burst12 /usr/local/bin/
COPY data /usr/local/src/phylogenize/
# RUN mkdir /usr/local/src/phylogenize
# COPY ./*.R /usr/local/src/phylogenize
WORKDIR /usr/local/src/phylogenize
RUN Rscript /usr/local/src/phylogenize/install-dependencies.R
CMD ["python", "aws-queue-reader.py"]
