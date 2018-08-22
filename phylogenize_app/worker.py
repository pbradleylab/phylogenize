#!/usr/bin/python3

import os, sys, shutil
import cgi
import uuid
import subprocess
from pystalkd.Beanstalkd import Connection
import re
import rpy2
import time

beanstalk = Connection(host='localhost', port=14711)
beanstalk.watch("phylogenize")
beanstalk.use("phylogenize")

MaxJobs = 2
JobList = [None] * MaxJobs
JobOutput = [None] * MaxJobs
JobErr = [None] * MaxJobs
JobTitle = [None] * MaxJobs
JobFile = [None] * MaxJobs
JobSlots = range(MaxJobs)

# Adapted from George V. Reilly:
# https://stackoverflow.com/questions/2032403/
def make_tarfile(output_filename, source_dir):
  files = os.listdir(source_dir)
  with tarfile.open(output_filename, "w:gz") as tar:
    for f in files:
      tar.add(os.path.join(source_dir, f),
          arcname=os.path.basename(source_dir))

while True:
  # poll for 5 seconds at a time, but don't bother polling if no open slots
  time.sleep(5)
  JobSlots = [i for i in range(MaxJobs) if JobList[i] == None]
  job = beanstalk.reserve(timeout = 0)
  if not job is None:
    print(("received job #%d out of %d: %s" % (MaxJobs - (len(JobSlots)) + 1, MaxJobs, job.body)))
    job_file_match = re.search("output_dir = \"([^\"]*)\"", job.body)
    if job_file_match is None:
      print("invalid job: %s" % job.body)
      job.delete()
      continue
    job_file = job_file_match.group(1)
    job_title = os.path.basename(os.path.dirname(job_file))
    if len(JobSlots) == 0:
      print("received job %s, but no open slots" % job_title)
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
        this_message = ("Waiting for the queue to be ready.\n\n" +
            "There are %d other jobs waiting." % (place_in_line))
        fh.write(this_message)
        fh.flush()
      print(this_message)
      print("Delaying %d seconds" % (jdelay))
      beanstalk.put(job.body, delay = jdelay)
    elif len(JobSlots) > 0:
      JobN = JobSlots[0]
      JobFile[JobN] = job_file
      JobOutput[JobN] = open(os.path.join(job_file, "progress.txt"), 'w')
      JobErr[JobN] = open(os.path.join(job_file, "stderr.txt"), 'w')
      JobList[JobN] = subprocess.Popen(["/usr/bin/Rscript", "-e", job.body],\
          stdout=JobOutput[JobN], \
          stderr=JobErr[JobN])
      JobTitle[JobN] = job_title
      JobOutput[JobN].flush()
      JobErr[JobN].flush()
      job.delete()
  for (N, sub) in enumerate(JobList):
    # Poll to see which have finished, then return them to None and clean up
    if not sub == None:
      JobOutput[N].flush()
      JobErr[N].flush()
      returncode = sub.poll()
      if not returncode == None:
        # job finished!
        print("Job %s finished running in slot %d - cleaning up" % \
            (JobTitle[N], N + 1))
        if os.path.isfile(os.path.join(job_file, "phylogenize-report.html")):
          make_tarfile(os.path.join(job_file, "output.tgz" % JobTitle[N]),
             os.path.join(job_file, "output"))
        JobErr[N].close()
        JobOutput[N].close()
        JobTitle[N] = None
        JobList[N] = None
        JobFile[N] = None
        print("...done!")
  # continue forever

