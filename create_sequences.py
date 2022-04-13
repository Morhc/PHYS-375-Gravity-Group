#controls everything

import os
import pandas as pd

from standards import Standards as s
from generate_star import create_star

def main():

    '''
    for num, Tc in enumerate(np.linspace(1e6, 50e6, 25), start=0):

        print(f'Generating star {num}')

        df, star_data = create_star(Tc)

        rho_c, Tc, R_star, M_star, L_star, T_star = star_data

        print('Saving the data')
        ### SAVING OUT THE DATA ###
        summary_file = os.path.join(s.data_folder, 'test', 'generated_stars.csv')

        with open(summary_file, 'a+') as f:
            f.write(f'{rho_c},{Tc},{R_star},{M_star},{L_star},{T_star},{num}\n')

        df.to_csv(rf"{os.path.join(s.data_folder, 'test', f'star_{num}.csv')}", index=False)

        del df
    '''

    Tc = 8235440
    df, star_data = create_star(Tc)

    rho_c, Tc, R_star, M_star, L_star, T_star = star_data

    summary_file = os.path.join(s.data_folder, 'summary_test.csv')

    with open(summary_file, 'a+') as f:
        f.write(f'{rho_c},{Tc},{R_star},{M_star},{L_star},{T_star}\n')

    df.to_csv(rf"{os.path.join(s.data_folder, 'test_star.csv')}", index=False)

    del df

if __name__ == "__main__":
    main()
