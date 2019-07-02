#!/usr/bin/env python3

import os, sys, shutil
import uuid
import subprocess
from pystalkd.Beanstalkd import Connection
import re
import time
import json
import tarfile
import click

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
          arcname=os.path.join('output/', f))

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
      print("received job %s, but no open slots" % jobdict["result_id"])
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
      with open(os.path.join(jobdict["output_dir"], "progress.txt"), "w") as fh:
        this_message = ("Waiting for the queue to be ready.\n\n" +
            "There are %d other jobs waiting." % (place_in_line))
        fh.write(this_message)
        fh.flush()
      print(this_message)
      print("Delaying %d seconds" % (jdelay))
      beanstalk.put(job.body, delay = jdelay)
      job.delete()
    elif len(JobSlots) > 0:
      # shutil.copy(os.path.join(jobdict["report_dir"], jobdict["report_name"]),
      #     jobdict["output_dir"])
      JobN = JobSlots[0]
      JobDict[JobN] = jobdict
      JobOutput[JobN] = open(os.path.join(jobdict["output_dir"], "progress.txt"), 'w')
      JobErr[JobN] = open(os.path.join(jobdict["output_dir"], "stderr.txt"), 'w')
      os.chdir(jobdict["output_dir"])
      # Don't use output_dir directly because then paths get screwed up on
      # render; also turn off cache which messes everything up
      Rcmd=(('sapply(c("phylogenize", "graphics", "stats", "methods",'
             '"grDevices", "biomformat"), function(.) library('
             'character.only=TRUE, .)); '
             'phylogenize::set_data_internal(); '
             'setwd("{output_dir}"); '
             'phylogenize::render.report('
             'do_cache=FALSE, '
             'output_file="{output_file}", '
             'input_format="{input_format}", '
             'biom_file="{biom_file}", '
             'db_version="{db_version}", '
             'abundance_file="{abundance_file}", '
             'metadata_file="{metadata_file}", '
             'separate_metadata={separate_metadata}, '
             'burst_cutoff="{burst_cutoff}", '
             'burst_dir="{burst_dir}", '
             'in_dir="{input_dir}", '
             'out_dir="{output_dir}", '
             'ncl={ncl}, '
             'type="{data_type}", '
             'which_phenotype="{which_phenotype}", '
             'which_envir="{which_envir}", '
             'dset_column="{dset_column}", '
             'env_column="{env_column}", '
             'sample_column="{sample_column}", '
             'assume_below_LOD={assume_below_lod_R}, '
             'single_dset={single_dset_R}, '
             'minimum={minimum}, '
             'relative_out_dir=".", '
             'working_dir="{output_dir}"'
             ')'
      ).format(
        output_file="index.html",
        input_dir=jobdict["input_dir"],
        output_dir=jobdict["output_dir"],
        burst_dir=jobdict["burst_dir"],
        ncl=1,
        separate_metadata=jobdict["separate_metadata"],
        db_version=jobdict["db_version"],
        abundance_file=jobdict["abundance_file"],
        metadata_file=jobdict["metadata_file"],
        biom_file=jobdict["biom_file"],
        data_type=jobdict["type"],
        input_format=jobdict["input_format"],
        which_phenotype=jobdict["which_phenotype"],
        which_envir=jobdict["which_envir"],
        dset_column=jobdict["dset_column"],
        env_column=jobdict["env_column"],
        sample_column=jobdict["sample_column"],
        burst_cutoff=format(jobdict["burst_cutoff"]),
        assume_below_lod_R="TRUE",
        single_dset_R=jobdict["single_dset"],
        minimum=jobdict["minimum"]
      ))
      job_rscript=Rcmd
      JobList[JobN] = subprocess.Popen(["/usr/bin/R", "-e", job_rscript], \
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
        if os.path.exists(JobDict[N]["input_dir"]) and JobDict[N]["cleanup"]:
          for filename in os.listdir(JobDict[N]["input_dir"]):
            rf = os.path.join(JobDict[N]["input_dir"], filename)
            print("Removing file %s..." % rf)
            os.remove(rf)
        if os.path.isfile(os.path.join(JobDict[N]["output_dir"], "index.html")):
          make_tarfile(os.path.join(JobDict[N]["output_dir"], "output.tgz"),
             JobDict[N]["output_dir"])
        JobErr[N].close()
        JobOutput[N].close()
        JobList[N] = None
        JobDict[N] = None
        print("...done!")
  # continue forever


