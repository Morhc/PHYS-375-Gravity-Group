#script: tools.py
#purpose: To have a hub for all useful functions.

import os
import pandas as pd
from standards import Standards as s

def load_data():
    """Loads in the data. Warning, the data is in CGMS units, NOT SI units.
    """

    highmass = pd.read_csv(os.path.join(s.data_folder, 'highmass_star.txt'), sep='\s+', header=0, skiprows=[1,2])
    lowmass = pd.read_csv(os.path.join(s.data_folder, 'lowmass_star.txt'), sep='\s+', header=0, skiprows=[1,2])
    summary = pd.read_csv(os.path.join(s.data_folder, 'star_summary.txt'), sep='\s+', header=0, skiprows=[1,2])

    return highmass, lowmass, summary
