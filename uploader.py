#!/usr/bin/python

import cgi
import cgitb; cgitb.enable()
import os, sys
import uuid
#from subprocess32 import run
import subprocess

print "Content-Type: text/html"
print
print """\
    <html>
    <body>
    <h2>Hello World!</h2>
    """

form = cgi.FieldStorage()
abd_f = form['abundance_file']
abd = abd_f.file.read()
metad_f = form['metadata_file']
mdf = metad_f.file.read()
dtype = form['dtype'].value
db = form['database'].value
phenotype = form['phenotype'].value
if phenotype == "provided":
  pheno_file = [l for l in form['phenotype_file'].file.readline()]
prior_type = form['prior_type']
if prior_type == "provided":
  prior_file = [l for l in form['phenotype_file'].file.readline()]
minimum = int(form['minimum'].value)
which_envir = form['which_envir'].value

### NOTE: got to sanitize all this

# make a new random directory

while True:
  direc = os.path.abspath(os.path.join(".", str(uuid.uuid4())))
  if not os.path.exists(direc):
    break

IDir = os.path.join(direc, "input")
ODir = os.path.join(direc, "output")
os.mkdir(direc)
os.mkdir(IDir)
os.mkdir(ODir)

with open(os.path.join(IDir, "abundance.tab"), 'w') as fh:
  for l in abd:
    fh.write(l)
with open(os.path.join(IDir, "metadata.tab"), 'w') as fh:
  for l in mdf:
    fh.write(l)

if phenotype == "provided":
  with open(os.path.join(IDir, "phenotype.tab"), 'w') as fh:
    for l in pheno_file:
      fh.write(l)
if prior_type == "provided":
  with open(os.path.join(IDir, "prior.tab"), 'w') as fh:
    for l in prior_file:
      fh.write(l)

rmark_render_cmd = "rmarkdown::render(\"plm-report-generator.Rmd\", " + \
    "output_format = \"html_document\", " + \
    ("output_dir = \"%s\", " % (ODir)) + \
    ("params = list(type = \"%s\", " % (dtype)) + \
    ("out_dir = \"%s\", " % (ODir)) + \
    ("in_dir = \"%s\", " % (IDir)) + \
    ("abundance_file = \"abundance.tab\", ") + \
    ("metadata_file = \"metadata.tab\", ") + \
    ("phenotype_file = \"phenotype.tab\", ") + \
    ("db_version = \"%s\", " % db) + \
    ("which_phenotype = \"%s\", " % phenotype) + \
    ("which_envir = \"%s\", " % which_envir) + \
    ("prior_type = \"%s\", " % prior_type) + \
    ("prior_file = \"prior.tab\", ") + \
    ("minimum = %d " % minimum) + \
    "))"

print(rmark_render_cmd)
print """
    </body>
    </html>
    """

os.environ["HOME"] = "/tmp/"
p = subprocess.Popen(["Rscript", "-e", rmark_render_cmd],\
    stdout=subprocess.PIPE, \
    stderr=subprocess.PIPE)
PipeOut = p.stdout
PipeErr = p.stderr
print(p.communicate())

