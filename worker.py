#!/usr/bin/python3

import os, sys, beanstalkc
import cgi
import uuid
import subprocess
from pystalkd.Beanstalkd import Connection
import re
import rpy2
import time

beanstalk = Connection(host='localhost', port=14711)
beanstalk.watch("phylogenize")

MaxJobs = 1
JobList = [None] * MaxJobs
JobOutput = [None] * MaxJobs
JobErr = [None] * MaxJobs
JobSlots = range(MaxJobs)

while True:
  # poll for 10 seconds at a time
  job = beanstalk.reserve(timeout = 10)
  # find out which slots are open
  JobSlots = [i for i in range(MaxJobs) if JobList[i] == None]
  if not job is None:
    print(("received job #%d out of %d: %s" % (JobN, MaxJobs, job.body)))
    job_file_match = re.search("output_dir = \"([^\"]*)\"", job.body)
    job_file = job_file_match.group(1)
    if len(JobSlots) > 0:
      JobN = JobSlots[0]
      JobOutput[JobN] = open(os.path.join(job_file, "progress.txt"), 'w')
      JobErr[JobN] = open(os.path.join(job_file, "stderr.txt"), 'w')
      JobList[JobN] = subprocess.Popen(["/usr/bin/Rscript", "-e", job.body],\
          stdout=JobOutput[JobN], \
          stderr=JobErr[JobN])
      job.delete()
    else:
      # no slots free, put it back on the stack and wait
      place_in_line = beanstalk.stats()["current-jobs-delayed"]
      if place_in_line > 10:
        jdelay = 600
      if place_in_line > 5:
        jdelay = 300
      if place_in_line > 2:
        jdelay = 150
      else:
        jdelay = 30
      with open(os.path.join(job_file, "progress.txt"), "w") as fh:
        fh.write("Waiting for the queue to be ready.\n\n" +
            "There are %d other jobs waiting." % place_in_line)
      beanstalk.put(job.body, delay = jdelay)
  for (N, sub) in enumerate(JobList):
    # Poll to see which have finished, then return them to None and clean up
    if not sub == None:
      returncode = sub.poll()
      if not returncode == None:
        # job finished!
        JobErr[N].close()
        JobOutput[N].close()
        JobList[N] = None


