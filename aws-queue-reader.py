import boto3
import os, sys
import cgi
import uuid
import subprocess
import re
import rpy2
import time

sqs = boto3.resource("sqs")
queue = sqs.get_queue_by_name(QueueName='phylogenize')

MaxJobs = 1
JobList = [None] * MaxJobs
JobOutput = [None] * MaxJobs
JobErr = [None] * MaxJobs
JobTitle = [None] * MaxJobs
JobSlots = range(MaxJobs)

while True:
  # poll for 5 seconds at a time, but don't bother polling if no open slots
  time.sleep(5)
  JobSlots = [i for i in range(MaxJobs) if JobList[i] == None]
  if len(JobSlots) > 0:
    # see if we can get a job going
    messages = queue.receive_messages(MaxNumberOfMessages=1)
    if len(messages) > 0:
      message = messages[0]
      print(("received job #%d out of %d: %s" % (
        MaxJobs - (len(JobSlots)) + 1,
        MaxJobs,
        message['Body'])))
      job_file_match = re.search("output_dir = \"([^\"]*)\"", message['Body'])
      job_file = job_file_match.group(1)
      job_title = os.path.basename(os.path.dirname(job_file))
      JobN = JobSlots[0]
      JobOutput[JobN] = open(os.path.join(job_file, "progress.txt"), 'w')
      JobErr[JobN] = open(os.path.join(job_file, "stderr.txt"), 'w')
      JobList[JobN] = subprocess.Popen(
        ["/usr/bin/Rscript",
          "-e",
          message['Body']],
        stdout=JobOutput[JobN],
        stderr=JobErr[JobN])
      JobTitle[JobN] = job_title
      JobOutput[JobN].flush()
      JobErr[JobN].flush()
      message.delete()
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
        JobErr[N].close()
        JobOutput[N].close()
        JobTitle[N] = None
        JobList[N] = None
        print("...done!")
    # continue forever

