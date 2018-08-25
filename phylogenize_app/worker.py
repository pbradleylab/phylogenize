#!/usr/bin/python3

import os, sys, shutil
import cgi
import uuid
import subprocess
from pystalkd.Beanstalkd import Connection
import re
import rpy2
import time
import json

beanstalk = Connection(host='localhost', port=14711)
beanstalk.watch("phylogenize")
beanstalk.use("phylogenize")

MaxJobs = 1
JobList = [None] * MaxJobs
JobOutput = [None] * MaxJobs
JobErr = [None] * MaxJobs
JobDict = [None] * MaxJobs
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
    # load json job
    try:
      jobdict = json.loads(job.body)
    except Exception as e:
      print("Exception %s: invalid job, %s" % (e, job.body))
      job.delete()
      continue
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
      with open(os.path.join(jobdict["ODir"], "progress.txt"), "w") as fh:
        this_message = ("Waiting for the queue to be ready.\n\n" +
            "There are %d other jobs waiting." % (place_in_line))
        fh.write(this_message)
        fh.flush()
      print(this_message)
      print("Delaying %d seconds" % (jdelay))
      beanstalk.put(job.body, delay = jdelay)
    elif len(JobSlots) > 0:
      shutil.copy(os.path.join(jobdict["report_dir"], jobdict["report_file"]),
          jobdict["ODir"])
      JobN = JobSlots[0]
      JobDict[JobN] = jobdict
      JobOutput[JobN] = open(os.path.join(jobdict["ODir"], "progress.txt"), 'w')
      JobErr[JobN] = open(os.path.join(jobdict["ODir"], "stderr.txt"), 'w')
      os.chdir(jobdict["ODir"])
      # Don't use output_dir directly because then paths get screwed up on
      # render
      job_rscript = \
        ("rmarkdown::render(\"%s\", " % (jobdict["report_file"]) +\
          "output_format = \"html_document\", " + \
          ("params = list(type = \"%s\", " % (jobdict["datatype"])) + \
          ("out_dir = \"%s\", " % (jobdict["ODir"])) + \
          ("in_dir = \"%s\", " % (jobdict["IDir"])) + \
          ("abundance_file = \"%s\", " % (jobdict["abundance_file"])) + \
          ("metadata_file = \"%s\", " % (jobdict["metadata_file"])) + \
          ("biom_file = \"%s\", " % (jobdict["biom_file"])) + \
          ("input_format = \"%s\", " % jobdict["input_format"]) + \
          ("phenotype_file = \"%s\", " % jobdict["phenotype_file"]) + \
          ("db_version = \"%s\", " % jobdict["database"]) + \
          ("which_phenotype = \"%s\", " % jobdict["phenotype"]) + \
          ("which_envir = \"%s\", " % jobdict["which_envir"]) + \
          ("prior_type = \"%s\", " % jobdict["prior_type"]) + \
          ("prior_file = \"%s\", " % jobdict["prior_file"]) + \
          ("minimum = %d" % int(jobdict["minimum"])) + \
          "))"
      JobList[JobN] = subprocess.Popen(["/usr/bin/Rscript", "-e", job_rscript],\
          stdout=JobOutput[JobN], \
          stderr=JobErr[JobN])
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
            (JobDict[N]["result_id"], N + 1))
        # clean up input files, which could be large and are no longer needed
        if os.path.exists(JobDict[N]["IDir"]):
          for filename in os.listdir(JobDict[N]["IDir"]):
            rf = os.path.join(JobDict[N]["IDir"], filename)
            print("Removing file %s..." % rf)
            os.remove(rf)
        if os.path.isfile(os.path.join(JobDict[N]["ODir"], "phylogenize-report.html")):
          make_tarfile(os.path.join(JobDict[N]["ODir"], "output.tgz"),
             JobDict[N]["ODir"])
        JobErr[N].close()
        JobOutput[N].close()
        JobList[N] = None
        JobDict[N] = None
        print("...done!")
  # continue forever


