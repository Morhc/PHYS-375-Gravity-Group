#controls everything

import os
import pandas as pd

from standards import Standards as s
from generate_star import create_star

def main():

    '''
    for num, scale in enumerate(range(1, 100+1), start=0):

        print(f'Generating star {num}')

        Tc = scale*1e6  # (core temperature in K) we will choose Tc to be the MS paramter

        df, star_data = create_star(Tc)

        rho_c, Tc, R_star, M_star, L_star, T_star = star_data

        print('Saving the data')
        ### SAVING OUT THE DATA ###
        summary_file = os.path.join(s.data_folder, 'stars_lam_10', 'generated_stars.csv')

        with open(summary_file, 'a+') as f:
            #converting to cgms units
            f.write(f'{rho_c},{Tc},{R_star},{M_star},{L_star},{T_star},{num}\n')

        df.to_csv(rf"{os.path.join(s.data_folder, 'stars_lam_10', f'star_{num}.csv')}", index=False)
    '''

    Tc = 8.23544e6
    df, star_data = create_star(Tc)
    rho_c, Tc, R_star, M_star, L_star, T_star = star_data
    summary_file = os.path.join(s.data_folder, 'summary_test.csv')
    with open(summary_file, 'a+') as f:
        #converting to cgms units
        f.write(f'{rho_c},{Tc},{R_star},{M_star},{L_star},{T_star}\n')
    df.to_csv(rf"{os.path.join(s.data_folder, 'test_star.csv')}", index=False)

if __name__ == "__main__":
    main()
