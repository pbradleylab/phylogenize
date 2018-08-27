#!./venv/bin/python3

NDAYS = 7

import os, sys, shutil, time, datetime
from datetime import date

curr_time = date.fromtimestamp(time.time())
last_week = curr_time - datetime.timedelta(days=NDAYS)

results_dir = os.environ.get('PHYLOGENIZE_RESULTS_DIR') or \
  "/home/pbradz/projects/phylogenize/instance/results/"

results_dir = os.path.abspath(results_dir)

results_contents = os.listdir(results_dir)

print("looking under %s" % (results_dir))

for rd in results_contents:
  p = os.path.join(results_dir, rd)
  if os.path.isdir(p):
    modified_recently = False
    for root, dirs, files in os.walk(p):
      for name in files:
        mtime = os.path.getmtime(os.path.join(root, name))
        dttime = date.fromtimestamp(mtime)
        if (dttime > last_week):
          modified_recently = True
          continue
      if modified_recently:
        continue
    if not modified_recently:
      if os.path.abspath(p).startswith(results_dir):
        print("removing: %s" % os.path.abspath(p))
        shutil.rmtree(p)

