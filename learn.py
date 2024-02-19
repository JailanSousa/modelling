import pandas as pd
import re

best = pd.read_csv('best_model.csv')

tst = [''.join(re.findall(r'.+(?=\.B)', x)) for x in best.name]
print(tst)
