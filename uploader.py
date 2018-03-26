#!/usr/bin/python3

import cgi
import cgitb; cgitb.enable()
import os, sys
import uuid
#from subprocess32 import run
import subprocess
from pystalkd.Beanstalkd import Connection
import bleach

beanstalk = Connection(host='localhost', port=14711)
beanstalk.use("phylogenize")

print("Content-Type: text/html")
print()
print("""\
    <html>
    <head><link rel="stylesheet" href="phylo.css"> </head>
    <body>
    <h2>Your job has been submitted</h2>
    """)
sys.stdout.flush()

# store and sanitize form output

form = cgi.FieldStorage()
abd_f = bleach.clean(form['abundance_file'])
abd = abd_f.file.read()
metad_f = bleach.clean(form['metadata_file'])
mdf = metad_f.file.read()
dtype = bleach.clean(form['dtype'].value)
db = form['database'].value
phenotype = bleach.clean(form['phenotype'].value)
if phenotype == "provided":
  pheno_file = [l for l in form['phenotype_file'].file.readline()]
prior_type = bleach.clean(form['prior_type'].value)
if prior_type == "provided":
  prior_file = [l for l in form['phenotype_file'].file.readline()]
minimum = int(bleach.clean(form['minimum'].value))
which_envir = bleach.clean(form['which_envir'].value)

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

rmark_render_cmd = "rmarkdown::render(\"phylogenize-report.Rmd\", " + \
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

#print(rmark_render_cmd)

print(("""
 Your output will appear <a href="./%s/output">here</a><p>
""" % (os.path.basename(direc))))
print(("""
<p>
<iframe width="50%%" src="%s" frameborder=2 id="progress"></iframe>
<script>
window.setInterval("reloadIFrame();", 30000);
function reloadIFrame() {
  document.getElementById("progress").src="%s"
}
</script>
<p>
    </body>
    </html>
    """ % (os.path.join(os.path.basename(direc), "output", "progress.txt"), \
    os.path.join(os.path.basename(direc), "output", "progress.txt"))))

sys.stdout.flush()

os.environ["HOME"] = "/tmp/"
beanstalk.put(rmark_render_cmd)

#p = subprocess.Popen(["Rscript", "-e", rmark_render_cmd],\
#    stdout=subprocess.PIPE, \
#    stderr=subprocess.PIPE)
#PipeOut = p.stdout
#PipeErr = p.stderr
#
#sys.stdout.flush()
#
#while True:
#  output = p.stdout.read(1)
#  if output == '' and p.poll() != None:
#    break
#  if output != '':
#    sys.stdout.write(output)
#    sys.stdout.write("\n")
#    sys.stdout.flush()
