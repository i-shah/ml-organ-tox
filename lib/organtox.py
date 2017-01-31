# 
#  Imran Shah 
#  imran.a.shah@gmail.com
#  

import pandas as pd
import numpy as np
import pymongo,re,math,pickle,copy,warnings
import sklearn.metrics as metrics
import ipyparallel as PP
from bson import json_util
from bson.objectid import ObjectId
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from matplotlib.colors import ListedColormap
from sklearn import neighbors, datasets
import pandas as pd
from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.datasets import make_classification
from sklearn import cross_validation
from sklearn.feature_selection import SelectKBest,f_classif
from sklearn.metrics import accuracy_score,f1_score,precision_score,recall_score
from sklearn.ensemble import RandomForestClassifier

warnings.simplefilter('ignore')


def flatten(L):
    if not L: return []
    if type(L) != list or len(L)==1: return L
    def aaa(a,b):
        if type(a)!=list: a=[a]
        if type(b)!=list: b=[b]
        return a+b

    return [i for i in list(set(reduce(aaa,L))) if i]

def initParallel(Data=None,Code=None,parallel_machine='your_parallel_machine'):
    RC = PP.Client(profile=parallel_machine)
    global lb_view
    global d_view
    d_view = RC[:]
    d_view.block = True
    lb_view = RC.load_balanced_view()
    lb_view.block = True
    
    if Code: 
        d_view.execute(Code)
    if Data:
        d_view.push(Data)

def openMongo(host="hostname",user='devel',passwd='devel',db='organtox_v1'):
    con2 = pymongo.MongoClient("mongodb://%s:%s@%s/%s" % (user,passwd,host,db))
    DB = con2[db]
    return DB


def getChemFP(CID,col=None,ds=None,fp=None,fill=None):
    Agg = [
            # Match chemicals in cluster
            {'$match': {
                     'dsstox_cid':{'$in':CID}}
            },
            # Include these fields
            {'$project':{'dsstox_cid':1,'_id':0,
                        'fp':'$'+ds},
            },
            # Unwind the fp 
            {'$unwind':"$fp"}
            ]
    if fp: Agg.append({'$match': {'fp':{'$in': fp}}})
    
    X = col.aggregate(Agg,allowDiskUse=True)
    if not X: return
    #try:
    #    R = pd.DataFrame(X['result'])
    #except:
    R = pd.DataFrame(list(X))

    if R.shape[0]==0 or R.shape[1]==0: return pd.DataFrame()
    R['x']=1
    return pd.pivot_table(R,index=['dsstox_cid'],columns='fp',values='x',
                          aggfunc=len,fill_value=fill)


def getToxDataSet(tox,MDB=None):
    if not MDB: return
    Pos = np.array([i['dsstox_cid'] for i in MDB.tox_fp.find({'tox_fpp1.ds':tox,'dsstox_cid':{'$exists':1}},{'_id':0,'dsstox_cid':1})])
    Neg = np.array([i['dsstox_cid'] for i in MDB.tox_fp.find({'tox_fpn1.ds':tox,'dsstox_cid':{'$exists':1}},{'_id':0,'dsstox_cid':1})])
    n_n = len(Neg)
    n_p = len(Pos)
    CID = list(np.concatenate((Pos,Neg)))
    X_tox = pd.DataFrame(np.concatenate((np.ones(n_p),np.zeros(n_n))),index=CID,columns=[tox])
    X_chm = getChemFP(CID,col=MDB.chm_fp,ds='mrgn.ds',fill=0)
    X_ct  = getChemFP(CID,col=MDB.chm_fp,ds='chmtp1.ds',fill=0)
    X_bio = getChemFP(CID,col=MDB.bio_fp,ds='bio1.ds',fill=0)
    X_bc  = pd.merge(X_bio,X_chm,left_index=True,right_index=True)
    X_bct = pd.merge(X_bio,X_ct,left_index=True,right_index=True)
    
    return dict(tox=X_tox,chm=X_chm,ct=X_ct,bio=X_bio,bc=X_bc,bct=X_bct)


def getDataSubSet(Y,X,n_ss1=None,n_ss2=None,missing='drop',max_xy=False):
    """
    n_ss1  : the number of positives / negatives for subset 1
    n_ss2 : the number of positives / negatives for subset 2
    """
    
    ID = Y.index.intersection(X.index)
    if ID.shape[0]==0: 
        print "Not enough overlap"
        return
    
    Y = Y.ix[ID]
    X = X.ix[ID]

    X0 = pd.merge(X,Y,left_index=True,right_index=True)
    
    if max_xy: return X0
    
    cls_ds=Y.columns[0]
    Neg = Y.index[Y[cls_ds]==0]
    Pos = Y.index[Y[cls_ds]==1]
    n_n = len(Neg)
    n_p = len(Pos)
  
    N    = n_p if n_p>n_n else n_n
    n_ss1 = n_ss1 if n_ss1<N else N

    # Subset 1
    I_p = np.random.randint(0,n_p,n_ss1)
    I_n = np.random.randint(0,n_n,n_ss1)
    CID1  = np.concatenate((Pos[I_p],Neg[I_n]))

    Cols=list(X.columns)+[cls_ds]
    X1 = X0.ix[CID1,Cols]
    
    # Subset 2
    P2 = list(set(range(n_p)).difference(I_p))
    N2 = list(set(range(n_n)).difference(I_n))
    
    if len(P2)>n_ss2 and len(N2)>n_ss2:
    
        I2_p = np.random.choice(P2,n_ss2,replace=False)
        I2_n = np.random.choice(N2,n_ss2,replace=False)

        CID2 = np.concatenate((Pos[I2_p],Neg[I2_n]))
        X2= X0.ix[CID2,Cols]

        return X1,X2
    
    else:
        
        return X1,pd.DataFrame()

def runOrganToxML(y_out,dt,n_ss1,r_seed,
                  n_ss2=0,n_ds_min=5,n_ds_max=100,n_ds_step=5,
                  ss_iters=10,cv_iters=10,cv_nfolds=10,
                  Col_ds=None,Col_lr=None,MDB=None
                  ):
    """
    run the ML strategy for y_out
    """
    np.random.seed(r_seed)
    DS = getToxDataSet(y_out,MDB)
    for j in range(ss_iters):
        Dss_cv,Dss_ev = getDataSubSet(DS['tox'][[y_out]],DS[dt],n_ss1,n_ss2)
        CLF_i = getClassifiers()
        for i_ds in range(n_ds_min,n_ds_max,n_ds_step):
            runML1(Dss_cv,y_out,CLF=CLF_i,
                   dt_out='tox',dt_in=dt,n_ds=i_ds,
                   cv_iters=cv_iters,cv_nfolds=cv_nfolds,
                   Col_ds=Col_ds,Col_lr=Col_lr)

    return True


def runML1(Data_cv,y,
           dt_out=None,dt_in=None,
           n_ds = 10,
           cv_nfolds=10,cv_iters=5,
           CLF=None,
           DS_top=None,
           save=False,
           Col_ds=None,
           Col_lr=None):
    """
    Data_cv: Data for cross validation
    y      : Column of Data_cv with classification attribute
    n_ds   : number of top descriptors to use
    """
    
    # Data For CV 
    X1,Y1 = Data_cv.drop(y,axis=1),Data_cv[y]
    ID_cv = X1.index
    cls_ds=y
    Neg = Data_cv.index[Data_cv[cls_ds]==0]
    Pos = Data_cv.index[Data_cv[cls_ds]==1]
    n_n = len(Neg)
    n_p = len(Pos)
    
    if not DS_top: DS_top=pd.DataFrame(np.zeros(X1.shape[1]),index=X1.columns,columns=['N'])

    
    # The dataset
    R_ds=dict(pred=y,dt_in=dt_in,dt_out=dt_out,n_ds=n_ds,n_obs=X1.shape[0],n_neg=n_n,n_pos=n_p,
              cvt=dict(pos=list(Pos),neg=list(Neg)))
    
    if Col_ds:
        ds_id = Col_ds.save(R_ds)
    
    Perf = []                        

    # Iterate through cross-validation
    for k_cv in range(cv_iters):
        P = runXVal(X1,Y1,n_ds=n_ds,cv_nfolds=cv_nfolds,DS_top=DS_top)
        Perf += P
        
    P1 =summarizePerf(Perf)
    P1['n_iter']=cv_iters
    P1a = {i['lr']:i for i in P1.to_dict('records')}
    
    # Store the full classifiers
    CLF=getClassifiers()
    FS   = SelectKBest(f_classif)  
    FS.set_params(k=n_ds)
    FS.fit_transform(X1,y=Y1)
    DS=X1.columns[FS.get_support()]
    DS_top.sort('N',ascending=False,inplace=True)
    DS=DS_top.ix[:n_ds].index

    LR_db={}
    for Nm,Clf in CLF.iteritems():
        Clf.fit(X1[DS],Y1)
        try:
            Y =Clf.predict(X1[DS])
        except:
            print " >> Training evaluation failed",Nm
        else:
            R=dict(pred=y,dt_in=dt_in,dt_out=dt_out,n_ds=n_ds,
                   lr=Nm,n_pos=n_p,n_neg=n_n,n_obs=X1.shape[0])
            P=dict(pt='all_obs')
            P.update(evalPred(Y1,Y))
            P['ds']=list(DS)
            R['perf_trn']=P
            R['perf_cvt']=P1a[Nm]
            #R['clf_pkl']=pickle.dumps(Clf)
            #R['qmrf']   = 'TBD'
            #R['ds_id'] = ds_id,
            LR_db[Nm]=Col_lr.insert(R)


def getClassifiers():
    Classifiers0 = dict(NB= GaussianNB(), 
                        KNN0=KNeighborsClassifier(3), 
                        KNN1=KNeighborsClassifier(algorithm='auto',n_neighbors=5,p=2,weights='uniform'), 
                        SVCL0=SVC(kernel='linear'), 
                        SVCR0=SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3,
                                  gamma=10.0, kernel='rbf', max_iter=100, probability=False, shrinking=True,
                                  tol=0.001, verbose=False),                    
                        CART0=DecisionTreeClassifier(max_depth=10),
                        CART1=DecisionTreeClassifier(max_features='auto'),
                        RF0=RandomForestClassifier()
                       )
    return Classifiers0

def evalPred(Y_truth,Y_inferred,post=None):
    M = dict(sens=recall_score(Y_truth,Y_inferred),
             spec=precision_score(Y_truth,Y_inferred),
             acc=accuracy_score(Y_truth,Y_inferred),
             f1=f1_score(Y_truth,Y_inferred))
       
    M['bacc']=0.5*(M['sens']+M['spec'])
   
    if post: M={k+post:v for k,v in M.iteritems()}
    
    return M

    
def summarizePerf(P):
    P0=pd.DataFrame(P)
    Agg = dict(bacc=dict(mn=np.mean,sd=np.std),
               f1=dict(mn=np.mean,sd=np.std),
               acc=dict(mn=np.mean,sd=np.std),
               sens=dict(mn=np.mean,sd=np.std),
               spec=dict(mn=np.mean,sd=np.std),
               n_train=dict(mn=np.median),
               n_test =dict(mn=np.median)
              )
    P1=P0.groupby(by=['lr','pt','cv_kfold']).aggregate(Agg)
    P2=P1.reset_index()
    C=['_'.join([k for k in i if k]) for i in P2.columns]
    P2.columns=C
    
    return P2

def runXVal(X1,Y1,n_ds=None,CLF=None,cv_nfolds=10,DS_top=pd.DataFrame()):
    Perf=[]
    FS   = SelectKBest(f_classif)    # Feature selection
    SKF = StratifiedKFold(Y1,n_folds=cv_nfolds,shuffle=True)
    if DS_top.shape[0]==0: DS_top=pd.DataFrame(np.zeros(X1.shape[1]),index=X1.columns)

    k_fold=0
    
    for I1,I2 in SKF:
        k_fold+=1
        I_train,I_test=X1.index[I1],X1.index[I2]
        FS.set_params(k=n_ds)
        FS.fit_transform(X1.ix[I_train],y=Y1.ix[I_train])

        DS=X1.columns[FS.get_support()]
        DS_top.ix[DS]+=1
            
        X_train,Y_train = X1.ix[I_train,DS],Y1.ix[I_train]
        X_test, Y_test  = X1.ix[I_test,DS], Y1.ix[I_test]
        CLF=getClassifiers()
        
        for Nm,Clf in CLF.iteritems():
            Clf.fit(X_train,Y_train)
            try:
                Y_pred=Clf.predict(X_test)
            except:
                print " >> Failed",Nm
            else:
                P=dict(lr=Nm,pt='cvt',n_train=len(I_train),n_test=len(I_test),
                       n_ds=n_ds,cv_kfold=cv_nfolds)
                P.update(evalPred(Y_test,Y_pred))
                Perf.append(P)

    return Perf

def predPerfSummary(tox,Col_ml=None,Col_sum=None,limit=1000):
    P_cvt=[]
    P_trn=[]

    for Pi in Col_ml.find({'pred':tox},
                          dict(_id=0,data=0,clf_pkl=0,qmrf=0)):
        x1 = Pi.pop('perf_trn')
        x2 = Pi.pop('perf_cvt')
        x1.pop('ds')
        x1.update(Pi)
        x2.update(Pi)

        P_trn.append(x1)
        P_cvt.append(x2)


    # CVT performance
    X1 = pd.DataFrame(P_cvt)
    agg1=dict(f1_mn=np.mean,bacc_mn=np.mean,acc_mn=np.mean,sens_mn=np.mean,spec_mn=np.mean)
    by0 = ['pred','lr','dt_in','n_obs','n_pos','n_neg','n_ds']
    Perf_cvt=X1.groupby(by=by0).aggregate(agg1)
    for x in ['acc','bacc', 'f1', 'sens', 'spec']:
        Perf_cvt['%s_sd' % x]=X1.groupby(by=by0).apply(lambda g: pooledStdev(g,x))
    
    Perf_cvt['pt']='cvt'
    Perf_cvt['dt_out']='tox'
    if Col_sum: Col_sum.insert_many(Perf_cvt.reset_index().to_dict('records'))
        
    # Training performance
    X1 = pd.DataFrame(P_trn)
    agg2={}
    for x in ['acc','bacc', 'f1', 'sens', 'spec']:
        agg2[x]={"%s_mn"%x:np.mean,"%s_sd"%x:np.std}

    Perf_trn=X1.groupby(by=by0).aggregate(agg2)
    Perf_trn.columns=[i[1] for i in Perf_trn.columns]
    Perf_trn = Perf_trn[sorted(Perf_trn.columns)]

    Perf_trn['dt_out']='tox'
    Perf_trn['pt']='trn'
    if Col_sum: Col_sum.insert_many(Perf_trn.reset_index().to_dict('records'))

def eta_squared(aov):
    aov['eta_sq'] = 'NaN'
    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])
    return aov
 
def omega_squared(aov):
    mse = aov['sum_sq'][-1]/aov['df'][-1]
    aov['omega_sq'] = 'NaN'
    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*mse))/(sum(aov['sum_sq'])+mse)
    return aov
 
