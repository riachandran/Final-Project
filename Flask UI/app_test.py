import os

from flask import Flask, request, redirect, render_template, jsonify, stream_with_context, Response
from werkzeug.utils import secure_filename
from multiprocessing import Pool, Process
from Bio import SeqIO

from pwm import pwm

app=Flask(__name__)

app.secret_key = "secret key"
app.config['MAX_CONTENT_LENGTH'] = 1000 * 1024 * 1024

path = os.getcwd()
# file Upload
UPLOAD_FOLDER = os.path.join(path, 'uploads')

if not os.path.isdir(UPLOAD_FOLDER):
    os.mkdir(UPLOAD_FOLDER)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

ALLOWED_EXTENSIONS = set(['fasta'])


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/')
def upload_form():
    return render_template('index.html')


@app.route('/uploadInputFile', methods=['POST','GET'])
def uploadInputFile():
    if request.method == 'POST':
        if 'file' not in request.files:
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            filename = 'filesubmited'
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            # return render_template('seqdisplay.html')
            return jsonify(loggedIn())
        else:
            return redirect(request.url) 
    elif request.method == "GET":
        return "ok" 
    
@app.route('/loggedIn', methods=['GET'])
def loggedIn():
    return render_template('seqdisplays.html')

@app.route('/processSeq', methods=['GET'])
def processSeq():
    def generateVectors():

        filename = 'filesubmited'
        input_file_name = os.path.join(UPLOAD_FOLDER, filename)
        sequences_dictionary = {sequence.id : sequence.seq for sequence in SeqIO.parse(input_file_name,'fasta')}
        attribute_data = [s.id for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
        seq_data = [s.seq for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
        final_feature_vector = [[] for _ in range(len(seq_data))]

        with Pool(os.cpu_count() - 2) as pool:
        # execute tasks in order
            for sequence,result in pool.imap_unordered(pwm, range(0,len(seq_data)),chunksize=30):
                final_feature_vector[sequence].extend(result)
                yield '<b>Proceesed sequence'+str(sequence)+':</b> '+str(seq_data[sequence]) +'<br/>\n'
    
    return Response(stream_with_context(generateVectors()), mimetype='text/html')


if __name__ == "__main__":
    app.run(host = '127.0.0.1',port = 5025, debug = False)