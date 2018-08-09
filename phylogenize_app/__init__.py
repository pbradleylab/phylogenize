import os
import shutil
import bleach
import uuid
import tempfile
import csv

from werkzeug.utils import secure_filename
from flask import Flask
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
    DataRequired
)
from wtforms import (
    FileField,
    RadioField,
    StringField
)
from wtforms.validators import Optional, InputRequired
from pystalkd.Beanstalkd import Connection
from flask_uploads import configure_uploads

RESULTS_FOLDER = './results/'
REPORT_TEMPLATE_PATH = './phylogenize-report.Rmd'
ALLOWED_EXTENSIONS = set([
  'txt',
  'tsv',
  'tab',
  'csv',
  'biom'
])
allowed_files = UploadSet('files', ALLOWED_EXTENSIONS)
beanstalk = Connection(host='localhost', port=14711)
beanstalk.use("phylogenize")

class JobForm(FlaskForm):
  abundances = FileField(validators = [
    Optional(),
    FileRequired(),
    FileAllowed(ALLOWED_EXTENSIONS, 'Only tab-delimited files or BIOM files accepted'
  ])
  metadata = FileField(validators = [
    Optional(),
    FileRequired(),
    FileAllowed(ALLOWED_EXTENSIONS, 'Only tab-delimited files or BIOM files accepted'
  ])
  biomfile = FileField(
    label="BIOM file",
    validators = [
      Optional(),
      FileRequired(),
      FileAllowed(ALLOWED_EXTENSIONS, 'Only tab-delimited files or BIOM files accepted'
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
    default = "midas_v1.0",
    label = "Database version"
  )
  phenotype = RadioField(
    choices = [
      ("prevalence", "prevalence"),
      ("specificity", "environmental specificity score")
    ]
  )
  which_envir = StringField("Environment", validators = InputRequired())

  def validate(self):
    if not super(JobForm, self).validate():
      return False
    if not (self.abundances.data and self.metadata.data) or biomfile.data:
      msg = 'Either two tab-delimited files (a species abundance table and a' +
      'sample annotation file) or a single BIOM file (including species' +
      'abundances and sample annotations) must be uploaded'
      self.abundances.errors.append(msg)
      self.metadata.errors.append(msg)
      self.biomfile.errors.append(msg)
      return False
    if self.abundances.data and self.metadata.data and biomfile.data:
      msg = 'Please only submit either two tab-delimited files (a species' +
      'abundance table and a sample annotation file) or a single BIOM file' +
      '(including species abundances and sample annotations), not both'
      self.abundances.errors.append(msg)
      self.metadata.errors.append(msg)
      self.biomfile.errors.append(msg)
      return False
    return True

def create_app(test_config=None):
  app = Flask(__name__, instance_relative_config=True)
  configure_uploads(app, allowed_files)
  app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY') or 'tankboozysurrealgrinning'
  app.config['WTF_CSRF_SECRET_KEY'] = os.environ.get('WTF_CSRF_SECRET_KEY') or\
    'tankboozysurrealgrinning'
  app.config['UPLOAD_FOLDER'] = os.environ.get('UPLOAD_FOLDER') or UPLOAD_FOLDER

  @app.route('/', methods=['GET', 'POST'])
  def home():
    form = JobForm()
    if request.method == 'POST':
      if form.validate_on_submit():
        new_result_id = process_form(form, request)
        return(redirect(url_for('results', result_id=new_result_id)))
    return(render_template('index.html', form = form))

  @app.route('/results/<result_id>')
  def display_results(result_id):
    result_id = secure_filename(result_id)
    direc = os.path.abspath(os.path.join(RESULTS_FOLDER, result_id))
    reportfile = os.path.join(direc,
        "output",
        "phylogenize_report.html")
    if os.file.exists(reportfile):
      completed = True
    else completed = False
    outfile = os.path.join(direc, "output", "progress.txt")
    errfile = os.path.join(direc, "output", "stderr.txt")
    return(render_template('results.html',
      completed = completed,
      result_id = result_id,
      outfile = outfile,
      errfile = errfile))

  @app.route('/results/<result_id>/<subfile>')
  def display_html(result_id):
    result_id = secure_filename(result_id)
    subfile = secure_filename(subfile)
    direc = os.path.abspath(os.path.join(RESULTS_FOLDER, result_id))
    if subfile == "phylogenize_report.html":
      reportfile = os.path.join(direc,
          "output",
          "phylogenize_report.html")
      if os.file.exists(reportfile):
        return(send_from_directory(direc, "phylogenize_report.html"))
      else:
        return(redirect(url_for('display_results'), result_id = result_id))
    elif subfile == "results.tgz":
      tgzfile = os.path.join(direc,
          "output",
          "results.tgz")
      if os.file.exists(reportfile):
        return(send_from_directory(direc, "results.tgz"))
      else:
        return(redirect(url_for('display_results'), result_id = result_id))
    else:
      return(redirect(url_for('display_results'), result_id = result_id))

  return(app)

def process_form(form = None, request = None):
  uploaded_biom = False
  if form.biomfile.data:
    uploaded_biom = True
  # Generate a unique UUID
  while True:
    new_result_id = str(uuid.uuid4())
    direc = os.path.abspath(os.path.join(RESULTS_FOLDER, new_result_id))
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
      os.path.join(IDir, "data.biom"))
  else:
    allowed_files.save(request.files['metadata'],
      os.path.join(IDir, "metadata.tab"))
    allowed_files.save(request.files['abundances'],
      os.path.join(IDir, "abundance.tab"))
  # Set some last options and submit to beanstalk queue
  minimum = 3
  prior_type = "uninformative"
  rmark_render_cmd = "rmarkdown::render(\"" + REPORT_TEMPLATE_PATH +
      "\", " +
      "output_format = \"html_document\", " + \
      ("output_dir = \"%s\", " % (ODir)) + \
      ("params = list(type = \"%s\", " % (dtype)) + \
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
