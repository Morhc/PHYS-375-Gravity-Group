#controls everything

import os
import pandas as pd

from standards import Standards as s
from generate_star import create_star

def main():
    Tc = 8.23544e+06  # (core temperature in K) we will choose Tc to be the MS paramter


    df, star_data = create_star(Tc)

    rho_c, Tc, R_star, M_star, L_star, T_star = star_data

    print('Saving the data')
    ### SAVING OUT THE DATA ###
    summary_file = os.path.join(s.data_folder, 'generated_stars.csv')

    with open(summary_file, 'a+') as f:
        #converting to cgms units
        f.write(f'{rho_c},{Tc},{R_star},{M_star},{L_star},{T_star}\n')

    df.to_csv(rf"{os.path.join(s.data_folder, 'star_1.csv')}", index=False)


if __name__ == "__main__":
    main()
