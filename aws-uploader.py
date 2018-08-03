#!/usr/bin/python3

import cgi
import cgitb; cgitb.enable()
import os, sys
import uuid
import bleach
import tempfile
import csv
import boto3

def sanitize_tsv(fhandle, outfile):
  reader = csv.reader(fhandle, delimiter = '\t')
  outfh = csv.writer(open(outfile, 'wb+'))
  writer = csv.writer(outfh, delimiter = '\t')
  header = reader.next()
  if (len(header) < 2):
    raise RuntimeError("Invalid tab-delimited file")
  header = [bleach.clean(x) for x in header]
  writer.writerow(header)
  for row in reader:
    rclean = [bleach.clean(x) for x in row]
    writer.writerow(rclean)
  outfh.close()

sqs = boto3.resource('sqs')
queue = sqs.create_queue(QueueName='phylogenize')

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
abd_f = form['abundance_file']
abd = abd_f.file
#abd = saferead(abd_f)
metad_f = form['metadata_file']
#mdf = saferead(metad_f)
mdf = metad_f.file
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
  direc = os.path.abspath(os.path.join("./results/", str(uuid.uuid4())))
  if not os.path.exists(direc):
    break

IDir = os.path.join(direc, "input")
ODir = os.path.join(direc, "output")
os.mkdir(direc)
os.mkdir(IDir)
os.mkdir(ODir)

#tsv_sanitize(abd_f.file, os.path.join(IDir, "abundance.tab"))
#tsv_sanitize(metad_f.file, os.path.join(IDir, "abundance.tab"))


with open(os.path.join(IDir, "abundance.tab"), 'wb') as fh:
  for l in abd:
    fh.write(l)
with open(os.path.join(IDir, "metadata.tab"), 'wb') as fh:
  for l in mdf:
    fh.write(l)

if phenotype == "provided":
  with open(os.path.join(IDir, "phenotype.tab"), 'wb') as fh:
    for l in pheno_file:
      fh.write(l)
if prior_type == "provided":
  with open(os.path.join(IDir, "prior.tab"), 'wb') as fh:
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
 Your output will appear <a href="./results/%s/output">here</a><p>
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
    """ % (os.path.join("results", os.path.basename(direc), "output", "progress.txt"), \
    os.path.join("results", os.path.basename(direc), "output", "progress.txt"))))

sys.stdout.flush()

os.environ["HOME"] = tempfile.mkdtemp()
queue.send_message(MessageBody=rmark_render_cmd)

