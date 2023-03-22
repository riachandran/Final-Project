import os
import timeit

from flask import Flask, request, redirect, render_template, jsonify, stream_with_context, Response
from werkzeug.utils import secure_filename
from multiprocessing import Pool, Process
from Bio import SeqIO
import numpy as np
import pandas as pd
from sklearn.model_selection import ShuffleSplit # or StratifiedShuffleSplit

from pwm import pwm
import models

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
        seq_data = [s.seq for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
        global final_feature_vector
        final_feature_vector = [[] for _ in range(len(seq_data))]

        with Pool(os.cpu_count() - 2) as pool:
        # execute tasks in order
            for sequence,result in pool.imap_unordered(pwm, range(0,len(seq_data)),chunksize=30):
                final_feature_vector[sequence].extend(result)
                yield '<b>Proceesed sequence '+str(sequence)+':</b> '+str(seq_data[sequence]) +'<br/>\n*******************************\n'

    return Response(stream_with_context(generateVectors()), mimetype='text/html')

@app.route('/runML', methods=['GET'])
def runML():
    return render_template('modelruns.html')

@app.route('/runRF', methods=['GET'])
def runRF():
    filename = 'filesubmited'
    input_file_name = os.path.join(UPLOAD_FOLDER, filename)
    attribute_data = [s.id for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
    host_names = attribute_data
    unique_hst = list(np.unique(host_names))
    
    unique_hst = [a.split('|')[6] for a in host_names]
    unique_hst = list(unique_hst)
    unique_virus = [a.split('|')[5] for a in host_names]
    unique_virus = list(np.unique(unique_virus))
    
    int_hosts = []
    for ind_unique in range(len(unique_hst)):
        variant_tmp = unique_hst[ind_unique]
        ind_tmp = unique_hst.index(variant_tmp)
        int_hosts.append(ind_tmp)

    print("Attribute data preprocessing Done")

    X = np.array(final_feature_vector)
    y = np.array(int_hosts)

    sss = ShuffleSplit(n_splits=1, test_size=0.3, random_state= 1)
    sss.get_n_splits(X, y)
    train_index, test_index = next(sss.split(X, y)) 

    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    start = timeit.default_timer()
    rf_return,prediction = models.rf_fun(X_train,y_train,X_test,y_test)
    stop = timeit.default_timer()
    print("RF Time : ", stop - start)
    keys = ["Accuracy","Precision","Recall","F1 (weighted)","F1 (Macro)","F1 (Micro)","ROC AUC"]

    rf_table_final = pd.DataFrame(
    {'Metrics': keys,
     'Result': rf_return
    })

    return render_template('rfresult.html', tables=[rf_table_final.to_html(index=False, table_id="tables1")], titles=['Metrics,Result'])

@app.route('/runSVM', methods=['GET'])
def runSVM():
    filename = 'filesubmited'
    input_file_name = os.path.join(UPLOAD_FOLDER, filename)
    attribute_data = [s.id for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
    host_names = attribute_data
    unique_hst = list(np.unique(host_names))
    
    unique_hst = [a.split('|')[6] for a in host_names]
    unique_hst = list(unique_hst)
    unique_virus = [a.split('|')[5] for a in host_names]
    unique_virus = list(np.unique(unique_virus))
    
    int_hosts = []
    for ind_unique in range(len(unique_hst)):
        variant_tmp = unique_hst[ind_unique]
        ind_tmp = unique_hst.index(variant_tmp)
        int_hosts.append(ind_tmp)

    print("Attribute data preprocessing Done")

    X = np.array(final_feature_vector)
    y = np.array(int_hosts)

    
    sss = ShuffleSplit(n_splits=1, test_size=0.3, random_state= 1)
    sss.get_n_splits(X, y)
    train_index, test_index = next(sss.split(X, y)) 

    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    start = timeit.default_timer()
    svm_return = models.svm_fun(X_train,y_train,X_test,y_test)
    stop = timeit.default_timer()
    print("SVM Time : ", stop - start)

    keys = ["Accuracy","Precision","Recall","F1 (weighted)","F1 (Macro)","F1 (Micro)","ROC AUC"]

    svm_table_final = pd.DataFrame(
    {'Metrics': keys,
     'Result': svm_return
    })

    return render_template('svmresult.html', tables=[svm_table_final.to_html(index=False, table_id="tables2")], titles=['Metrics,Result'])

@app.route('/runMLP', methods=['GET'])
def runMLP():
    filename = 'filesubmited'
    input_file_name = os.path.join(UPLOAD_FOLDER, filename)
    attribute_data = [s.id for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
    host_names = attribute_data
    unique_hst = list(np.unique(host_names))
    
    unique_hst = [a.split('|')[6] for a in host_names]
    unique_hst = list(unique_hst)
    unique_virus = [a.split('|')[5] for a in host_names]
    unique_virus = list(np.unique(unique_virus))
    
    int_hosts = []
    for ind_unique in range(len(unique_hst)):
        variant_tmp = unique_hst[ind_unique]
        ind_tmp = unique_hst.index(variant_tmp)
        int_hosts.append(ind_tmp)

    print("Attribute data preprocessing Done")

    X = np.array(final_feature_vector)
    y = np.array(int_hosts)

    
    sss = ShuffleSplit(n_splits=1, test_size=0.3, random_state= 1)
    sss.get_n_splits(X, y)
    train_index, test_index = next(sss.split(X, y)) 

    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    start = timeit.default_timer()
    mlp_return = models.mlp_fun(X_train,y_train,X_test,y_test)
    stop = timeit.default_timer()
    print("MLP Time : ", stop - start)

    keys = ["Accuracy","Precision","Recall","F1 (weighted)","F1 (Macro)","F1 (Micro)","ROC AUC"]

    mlp_table_final = pd.DataFrame(
    {'Metrics': keys,
     'Result': mlp_return
    })

    return render_template('mlpresult.html', tables=[mlp_table_final.to_html(index=False, table_id="tables3")], titles=['Metrics,Result'])

@app.route('/runKNN', methods=['GET'])
def runKNN():
    filename = 'filesubmited'
    input_file_name = os.path.join(UPLOAD_FOLDER, filename)
    attribute_data = [s.id for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
    host_names = attribute_data
    unique_hst = list(np.unique(host_names))
    
    unique_hst = [a.split('|')[6] for a in host_names]
    unique_hst = list(unique_hst)
    unique_virus = [a.split('|')[5] for a in host_names]
    unique_virus = list(np.unique(unique_virus))
    
    int_hosts = []
    for ind_unique in range(len(unique_hst)):
        variant_tmp = unique_hst[ind_unique]
        ind_tmp = unique_hst.index(variant_tmp)
        int_hosts.append(ind_tmp)

    print("Attribute data preprocessing Done")

    X = np.array(final_feature_vector)
    y = np.array(int_hosts)

    
    sss = ShuffleSplit(n_splits=1, test_size=0.3, random_state= 1)
    sss.get_n_splits(X, y)
    train_index, test_index = next(sss.split(X, y)) 

    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    start = timeit.default_timer()
    knn_return = models.knn_fun(X_train,y_train,X_test,y_test)
    stop = timeit.default_timer()
    print("KNN Time : ", stop - start)

    keys = ["Accuracy","Precision","Recall","F1 (weighted)","F1 (Macro)","F1 (Micro)","ROC AUC"]

    knn_table_final = pd.DataFrame(
    {'Metrics': keys,
     'Result': knn_return
    })

    return render_template('knnresult.html', tables=[knn_table_final.to_html(index=False, table_id="tables4")], titles=['Metrics,Result'])

@app.route('/runLR', methods=['GET'])
def runLR():
    filename = 'filesubmited'
    input_file_name = os.path.join(UPLOAD_FOLDER, filename)
    attribute_data = [s.id for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
    host_names = attribute_data
    unique_hst = list(np.unique(host_names))
    
    unique_hst = [a.split('|')[6] for a in host_names]
    unique_hst = list(unique_hst)
    unique_virus = [a.split('|')[5] for a in host_names]
    unique_virus = list(np.unique(unique_virus))
    
    int_hosts = []
    for ind_unique in range(len(unique_hst)):
        variant_tmp = unique_hst[ind_unique]
        ind_tmp = unique_hst.index(variant_tmp)
        int_hosts.append(ind_tmp)

    print("Attribute data preprocessing Done")

    X = np.array(final_feature_vector)
    y = np.array(int_hosts)

    
    sss = ShuffleSplit(n_splits=1, test_size=0.3, random_state= 1)
    sss.get_n_splits(X, y)
    train_index, test_index = next(sss.split(X, y)) 

    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    start = timeit.default_timer()
    lr_return = models.lr_fun(X_train,y_train,X_test,y_test)
    stop = timeit.default_timer()
    print("LR Time : ", stop - start)

    keys = ["Accuracy","Precision","Recall","F1 (weighted)","F1 (Macro)","F1 (Micro)","ROC AUC"]

    lr_table_final = pd.DataFrame(
    {'Metrics': keys,
     'Result': lr_return
    })

    return render_template('lrresult.html', tables=[lr_table_final.to_html(index=False, table_id="tables5")], titles=['Metrics,Result'])

@app.route('/runDT', methods=['GET'])
def runDT():
    filename = 'filesubmited'
    input_file_name = os.path.join(UPLOAD_FOLDER, filename)
    attribute_data = [s.id for s in SeqIO.parse(input_file_name,'fasta') if len(s.id.split('|')) > 6]
    host_names = attribute_data
    unique_hst = list(np.unique(host_names))
    
    unique_hst = [a.split('|')[6] for a in host_names]
    unique_hst = list(unique_hst)
    unique_virus = [a.split('|')[5] for a in host_names]
    unique_virus = list(np.unique(unique_virus))
    
    int_hosts = []
    for ind_unique in range(len(unique_hst)):
        variant_tmp = unique_hst[ind_unique]
        ind_tmp = unique_hst.index(variant_tmp)
        int_hosts.append(ind_tmp)

    print("Attribute data preprocessing Done")

    X = np.array(final_feature_vector)
    y = np.array(int_hosts)

    
    sss = ShuffleSplit(n_splits=1, test_size=0.3, random_state= 1)
    sss.get_n_splits(X, y)
    train_index, test_index = next(sss.split(X, y)) 

    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    start = timeit.default_timer()
    dt_return = models.fun_decision_tree(X_train,y_train,X_test,y_test)
    stop = timeit.default_timer()
    print("DT Time : ", stop - start)

    keys = ["Accuracy","Precision","Recall","F1 (weighted)","F1 (Macro)","F1 (Micro)","ROC AUC"]

    dt_table_final = pd.DataFrame(
    {'Metrics': keys,
     'Result': dt_return
    })

    return render_template('dtresult.html', tables=[dt_table_final.to_html(index=False, table_id="tables6")], titles=['Metrics,Result'])


if __name__ == "__main__":
    app.run(host = '127.0.0.1',port = 5109, debug = False)