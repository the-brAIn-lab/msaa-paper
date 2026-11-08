"""Normalised scans (intact, word, rest) as X.npz for fig5_per_seed.py. Writes outputs/multiseed_stability/X.npz."""
import numpy as np, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402
ROOT = paths.REPO
from scipy.io import loadmat
sys.path.insert(0, str(ROOT / 'src'))
from msaa_pipeline.data import _normalize_across_fallback
p=loadmat(paths.DATA)
X=[]; lab=[]
for c in ["intact","word","rest"]:
    arr=p[c]; subs=[arr[:,i][0] for i in range(arr.shape[1])]
    if c=="rest": subs=[np.asarray(x)[:-100] for x in subs]
    n=_normalize_across_fallback(subs); X+=n; lab+=[c]*len(n)
paths.X_NPZ.parent.mkdir(parents=True, exist_ok=True)
np.savez(paths.X_NPZ, *X, lab=np.array(lab))
print(len(X), X[0].shape, X[-1].shape, '->', paths.X_NPZ)
