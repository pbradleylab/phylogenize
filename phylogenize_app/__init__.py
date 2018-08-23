import os
import shutil
import bleach
import uuid
import tempfile
import csv
import re
import tarfile

from functools import wraps, update_wrapper
from werkzeug.utils import secure_filename
from datetime import datetime
from flask import Flask, make_response, send_from_directory
from flask import (
        Blueprint,
        flash,
        g,
        redirect,
        render_template,
        request,
        session,
        url_for
)
from flask_wtf import FlaskForm
from flask_wtf.file import (
    FileRequired,
    DataRequired,
    FileAllowed
)
from flask_wtf.recaptcha import RecaptchaField

from wtforms import (
    FileField,
    RadioField,
    StringField
)
from wtforms.validators import Optional, InputRequired
from pystalkd.Beanstalkd import Connection
from flask_uploads import configure_uploads, UploadSet

#RECAPTCHA_PUBLIC_KEY = '6LecpGkUAAAAAIy-lVILcUA-f6IsElHD-nn9OGuY'
RECAPTCHA_PUBLIC_KEY = '6LeIxAcTAAAAAJcZVRqyHh71UMIEGNQ_MXjiZKhI'
#RECAPTCHA_PRIVATE_KEY = '6LecpGkUAAAAAFuTGxJ6Um3wztmZxXASZBKKOwHi'
RECAPTCHA_PRIVATE_KEY = '6LeIxAcTAAAAAGG-vFI1TnRWxMZNFuojJ4WifJWe'
UPLOAD_FOLDER = os.path.abspath('/var/www/phylogenize/instance/results/')
REPORT_TEMPLATE_PATH = os.path.abspath('/var/www/phylogenize/phylogenize-report.Rmd')
ALLOWED_EXTENSIONS = set([
  'txt',
  'tsv',
  'tab',
  'csv',
  'biom'
])
beanstalk = Connection(host='localhost', port=14711)
beanstalk.use("phylogenize")

class JobForm(FlaskForm):
  abundances = FileField(validators = [
    Optional(),
    FileRequired(),
    FileAllowed(ALLOWED_EXTENSIONS,
      message='Only tab-delimited files or BIOM files accepted')
  ])
  metadata = FileField(validators = [
    Optional(),
    FileRequired(),
    FileAllowed(ALLOWED_EXTENSIONS,
      message='Only tab-delimited files or BIOM files accepted')
  ])
  biomfile = FileField(
    label="BIOM file",
    validators = [
      Optional(),
      FileRequired(),
      FileAllowed(ALLOWED_EXTENSIONS,
        message='Only tab-delimited files or BIOM files accepted')
    ]
  )
  datatype = RadioField(
    choices = [
      ("16S", "16S"),
      ("midas", "Shotgun (MIDAS)")
    ],
    default = "16S", 
    label = "Data type"
  )
  database = RadioField(
    choices = [
      ("midas_v1.0", "MIDAS 1.0"),
      ("midas_v1.2", "MIDAS 1.2")
    ],
    default = "midas_v1.2",
    label = "Database version"
  )
  phenotype = RadioField(
    choices = [
      ("prevalence", "prevalence"),
      ("specificity", "specificity")
    ], 
    default = "prevalence"
  )
  which_envir = StringField("Environment",
      validators = [InputRequired(message = 'Must provide an environment')],
      render_kw = {"placeholder": "e.g., stool"})

  def validate(self):
    if not super(JobForm, self).validate():
      return False
    if not ((self.abundances.data and self.metadata.data) or \
        self.biomfile.data):
      msg = "Either two tab-delimited files (a species abundance table and a" +\
        " sample annotation file) or a single BIOM file (including species" +\
        " abundances and sample annotations) must be uploaded"
      self.abundances.errors.append(msg)
      return False
    if self.abundances.data and self.metadata.data and self.biomfile.data:
      msg = "Please only submit either two tab-delimited files (a species" +\
        " abundance table and a sample annotation file) or a single BIOM file" +\
        " (including species abundances and sample annotations), not both"
      self.abundances.errors.append(msg)
      return False
    return True

def nocache(view):
  @wraps(view)
  def no_cache(*args, **kwargs):
    response = make_response(view(*args, **kwargs))
    response.headers['Last-Modified'] = datetime.now()
    response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, post-check=0, pre-check=0, max-age=0'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '-1'
    return response
  return update_wrapper(no_cache, view)

def create_app(config=None):

  app = Flask(__name__)
  app.config.from_pyfile('../phylogenize_default.cfg')
  app.config.from_pyfile('../instance/phylogenize.cfg', silent=True)

  allowed_files = UploadSet('files', app.config['ALLOWED_EXTENSIONS'])
  configure_uploads(app, allowed_files)

  @app.route('/', methods=['GET', 'POST'])
  @nocache
  def home():
    form = JobForm()
    if request.method == 'POST':
      if form.validate():
        new_result_id = process_form(allowed_files,
            form,
            request,
            os.path.join(app.config['APPLICATION_ROOT'],
              app.config['UPLOAD_FOLDER']))
        return(redirect(url_for('display_results', result_id=new_result_id)))
      flash("Error with form submission")
      return(render_template('index.html', form=form))
    return(render_template('index.html', form=form))

  @app.route('/tutorial', methods=['GET', 'POST'])
  @nocache
  def tutorial():
    return(render_template('tutorial.html'))

  @app.route('/about', methods=['GET', 'POST'])
  @nocache
  def about():
    return(render_template('about.html'))

  @app.route('/results')
  def reroute():
    return(redirect(url_for('home')))

  @app.route('/results/<result_id>')
  def display_results(result_id):
    result_id = secure_filename(result_id)
    direc = os.path.abspath(os.path.join(app.config['APPLICATION_ROOT'], app.config['UPLOAD_FOLDER'], result_id))
    reportfile = os.path.join(direc,
        "output",
        "phylogenize-report.html")
    if os.path.isfile(reportfile):
      completed = True
    else:
      completed = False
    outfile = os.path.join(direc, "output", "progress.txt")
    errfile = os.path.join(direc, "output", "stderr.txt")
    errmsgfile = os.path.join(direc, "output", "errmsg.txt")
    phylo_errortext = ''
    if os.path.isfile(errfile):
      with open(errfile, 'r') as fh:
        errlines = [l for l in fh.readlines()]
        if re.search(r'Execution halted', errlines[-1]):
          errormsg = True
          if os.path.isfile(errmsgfile):
            with open(errmsgfile, 'r') as fh:
              phylo_errortext = str("".join([l for l in fh.readlines()))
        else:
          errormsg = False
        if (len(errlines) >= 20) and not errormsg:
          errtext = str("".join(errlines[-20:]))
        else:
          errtext = str("".join(errlines))
    else:
      errtext = ''
      errormsg = False
    if os.path.isfile(outfile):
      with open(outfile, 'r') as fh:
        outlines = [l for l in fh.readlines()]
        if (len(outlines) >= 20) and not errormsg:
          outtext = str("".join(outlines[-20:]))
        else:
          outtext = str("".join(outlines))
        pctlines = list(filter(lambda x: re.search(r'\\|.*\\| *\%', x),
          outlines))
        pctline = pctlines[-1]
        pct = float(re.sub('.*\\| (.*)%', '\\1', pctline))
    else:
      outtext = ''
      pct = float(0)
    return(render_template('results.html',
      completed = completed,
      errormsg = errormsg,
      result_id = result_id,
      out = outtext,
      err = errtext,
      phylo_err = phylo_errortext,
      pct = pct,
      reportfile = reportfile))

  @app.route('/results/<result_id>/<subfile>')
  def display_result_file(result_id, subfile):
    result_id = secure_filename(result_id)
    subfile = secure_filename(subfile)
    direc = os.path.abspath(os.path.join(app.config['APPLICATION_ROOT'], app.config['UPLOAD_FOLDER'], result_id))
    if subfile == "phylogenize-report.html":
      reportfile = os.path.join(direc,
          "output",
          "phylogenize-report.html")
      if os.path.isfile(reportfile):
        return(send_from_directory(os.path.join(direc, "output"), "phylogenize-report.html"))
      else:
        return(redirect(url_for('display_results', result_id=result_id)))
    elif subfile == "output.tgz":
      tgzfile = os.path.join(direc,
          "output",
          "output.tgz")
      if os.path.isfile(tgzfile):
        return(send_from_directory(os.path.join(direc, "output"), "output.tgz"))
      else:
        return(redirect(url_for('display_results', result_id=result_id)))
    else:
      return(redirect(url_for('display_results', result_id=result_id)))

  return(app)

def process_form(allowed_files, form=None, request=None, upload_folder="../instance/results"):

  uploaded_biom = False
  if form.biomfile.data:
    uploaded_biom = True

  # Generate a unique UUID
  while True:
    new_result_id = str(uuid.uuid4())
    direc = os.path.abspath(os.path.join(upload_folder, new_result_id))
    if not os.path.exists(direc):
      break

  # Set up directory structure
  IDir = os.path.join(direc, "input")
  ODir = os.path.join(direc, "output")
  os.mkdir(direc)
  os.mkdir(IDir)
  os.mkdir(ODir)

  # Bleach form components
  database = str(bleach.clean(form.database.data))
  datatype = str(bleach.clean(form.datatype.data))
  phenotype = str(bleach.clean(form.phenotype.data))
  which_envir = str(bleach.clean(form.which_envir.data))
  input_format = ("biom" if uploaded_biom else "tabular")

  # Upload the BIOM or tabular files
  if uploaded_biom:
    allowed_files.save(request.files['biomfile'],
      folder = IDir,
      name = "data.biom")
  else:
    allowed_files.save(request.files['metadata'],
      folder = IDir,
      name = "metadata.tab")
    allowed_files.save(request.files['abundances'],
      folder = IDir,
      name = "abundance.tab")

  # Set some last options and submit to beanstalk queue
  minimum = 3
  prior_type = "uninformative"
  rmark_render_cmd = \
    ("rmarkdown::render(\"%s\", " % (REPORT_TEMPLATE_PATH)) +\
      "output_format = \"html_document\", " + \
      ("output_dir = \"%s\", " % (ODir)) + \
      ("params = list(type = \"%s\", " % (datatype)) + \
      ("out_dir = \"%s\", " % (ODir)) + \
      ("in_dir = \"%s\", " % (IDir)) + \
      ("abundance_file = \"abundance.tab\", ") + \
      ("metadata_file = \"metadata.tab\", ") + \
      ("biom_file = \"data.biom\", ") + \
      ("input_format = \"%s\", " % input_format) + \
      ("phenotype_file = \"phenotype.tab\", ") + \
      ("db_version = \"%s\", " % database) + \
      ("which_phenotype = \"%s\", " % phenotype) + \
      ("which_envir = \"%s\", " % which_envir) + \
      ("prior_type = \"%s\", " % prior_type) + \
      ("prior_file = \"prior.tab\", ") + \
      ("minimum = %d " % minimum) + \
      "))"
  beanstalk.put(rmark_render_cmd)
  return(new_result_id)
