#!/usr/bin/python3

import os, sys, beanstalkc
import cgi
import uuid
#from subprocess32 import run
import subprocess
from pystalkd.Beanstalkd import Connection
import re

beanstalk = Connection(host='localhost', port=14711)
beanstalk.watch("phylogenize")

while True:
  job = beanstalk.reserve()
  print(("received job %s" % job.body))
  job_file_match = re.search("output_dir = \"([^\"]*)\"", job.body)
  job_file = job_file_match.group(1)
  with open(os.path.join(job_file, "progress.txt"), 'w') as fh:
    with open(os.path.join(job_file, "stderr.txt"), 'w') as fh2:
      p = subprocess.Popen(["Rscript", "-e", job.body],\
          stdout=fh, \
          stderr=fh2)
  job.delete()
