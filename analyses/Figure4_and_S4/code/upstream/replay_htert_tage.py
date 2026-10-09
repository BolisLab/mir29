"""Optional exact tAge replay from frozen YuGene/control-centered expression.

Run from the bundle root with scikit-learn 1.3.2, pandas, numpy and joblib.
No model training, network requests or dependency installation occurs here.
"""
from pathlib import Path
import sys
import hashlib
import numpy as np
import pandas as pd
import sklearn

root = Path.cwd()
sys.path.insert(0, str(root / 'code/upstream/tage'))
import tage_predict

if sklearn.__version__ != '1.3.2':
    raise RuntimeError('The historical prediction environment requires scikit-learn 1.3.2.')
model_path = root / 'data/models/tAge_YuGene.pkl'
assert hashlib.md5(model_path.read_bytes()).hexdigest() == 'ad893ce7c008094b9f36450257d1934a'
matrix = pd.read_csv(root / 'data/models/htert_yugene_centered.tsv', sep='\t', index_col=0)
matrix.index = matrix.index.map(str)
model, features = tage_predict._load_clock(model_path)
aligned = tage_predict._align_features(matrix.T, features)
predicted = model.predict(aligned)
assert np.isfinite(predicted).all()
expected = pd.read_csv(root / 'data/models/htert_tAge_expected.tsv', sep='\t').set_index('sample')
np.testing.assert_allclose(predicted, expected.loc[matrix.columns, 'tAge_raw'], rtol=1e-10, atol=1e-12)
pd.DataFrame({'sample':matrix.columns,'tAge_raw':predicted}).to_csv(
    root / 'validation/htert_tAge_replayed.tsv', sep='\t', index=False)
print('PASS: all 42 hTERT tAge predictions match the frozen YuGene analysis.')
