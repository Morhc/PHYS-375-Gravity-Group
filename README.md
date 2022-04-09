# PHYS-375-Gravity-Group
This is the repository for the PHYS 375 group project where we will be modifying gravity.

This group will modify gravity in a specific way and infer the limits that may be placed on the modification length scale given that the Main Sequence is known to better than 10%. While many potential modifications of gravity may be considered, we will consider only one simple possible change: Below some scale λ we imagine that the gravitational force has a component that scales as r^-3.
The goal is then to constrain λ. Note that λ can be both positive and negative, and in this case will produce a linear cusp of one sort or another at the center of the star. The group should generate a variety of main sequences as functions of λ and explore the implications for the HR diagram, and the L-M and M-R relations.

## Dependencies
This project will require the user to install the following Python libraries:
- numpy: https://pypi.org/project/numpy/
- pandas: https://pypi.org/project/pandas/
- matplotlib: https://pypi.org/project/matplotlib/

**You may run into some issues with LaTeX.**

## How to Run
1. Set the lambda value you want in FofLambda
2. Set the Tc you want in create_sequences.py
3. Run create_sequences.py
