#controls everything

import os
import pandas as pd
import numpy as np

import time

from standards import Standards as s
from generate_star import create_star

lam = '-1e3'

def main():

    for num, Tc in enumerate(np.linspace(3e6, 25e6, 25)[0:], start=0):
        start = time.time()

        print(f'Generating star {num}. Temperature of {Tc}.')

        df, star_data = create_star(Tc)

        rho_c, Tc, R_star, M_star, L_star, T_star = star_data

        print(f'Saving the data. Took {time.time()-start} seconds.')
        ### SAVING OUT THE DATA ###
        summary_file = os.path.join(s.data_folder, f'stars_lam_{lam}', 'generated_stars.csv')

        with open(summary_file, 'a+') as f:
            f.write(f'{rho_c},{Tc},{R_star},{M_star},{L_star},{T_star},{num}\n')

        df.to_csv(rf"{os.path.join(s.data_folder, f'stars_lam_{lam}', f'star_{num}.csv')}", index=False)

        del df

    '''
    Tc = 8235440
    Tc = 2e7
    Tc = 3.5e7
    df, star_data = create_star(Tc)

    rho_c, Tc, R_star, M_star, L_star, T_star = star_data

    summary_file = os.path.join(s.data_folder, 'summary_test.csv')

    with open(summary_file, 'a+') as f:
        f.write(f'{rho_c},{Tc},{R_star},{M_star},{L_star},{T_star}\n')

    df.to_csv(rf"{os.path.join(s.data_folder, 'test_star.csv')}", index=False)

    del df
    '''

if __name__ == "__main__":
    main()
