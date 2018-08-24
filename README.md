# phylogenize


## Running locally without the web server

To run phylogenize locally without the web server, you can simply run the RMarkdown notebook "phylogenize-report.Rmd" and change or override the values in the YAML header. Look at 

beanstalkd -l 127.0.0.1 -p 14711
sudo -u www-data python3 worker.py

