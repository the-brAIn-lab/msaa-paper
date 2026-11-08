import numpy as np, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
from scipy.io import loadmat
sys.path.insert(0, str(ROOT / 'src'))
from msaa_pipeline.data import _normalize_across_fallback
p=loadmat(ROOT / 'data/input/pieman_data.mat')
X=[]; lab=[]
for c in ["intact","word","rest"]:
    arr=p[c]; subs=[arr[:,i][0] for i in range(arr.shape[1])]
    if c=="rest": subs=[np.asarray(x)[:-100] for x in subs]
    n=_normalize_across_fallback(subs); X+=n; lab+=[c]*len(n)
np.savez('X.npz', *X, lab=np.array(lab))
print(len(X), X[0].shape, X[-1].shape)
