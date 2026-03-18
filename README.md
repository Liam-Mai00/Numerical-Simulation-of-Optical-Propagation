# Numerical Simulation of Optical Propagation Python Repository
The repository aims to **simulate optical propagation in turbulent media**.\
This is based on Jason Schmidt's book on 'Numerical Simulation of Optical Propagation'. The book can be found on the [SPIE digital library](https://www.spiedigitallibrary.org/ebooks/PM/Numerical-Simulation-of-Optical-Wave-Propagation-with-Examples-in-MATLAB/eISBN-9780819483270/10.1117/3.866274)
## Repository Purpose
This repo is simply for my own eduction into Github, Python, and Optics. I am not affiliated with J. Schmidt but I am thankful of his work!
## Setup
The initial setup is recommended to prevent relative import errors and avoids
```python
sys.path.insert()
```
hacks.
[See here](https://stackoverflow.com/questions/6323860/sibling-package-imports/50193944#50193944)
- Open terminal
- Activate virtual environment
- Navigate to the repository
- Run
```python
pip install -e .
```
(editable format)
- Should see opt-prop package in site-package directory of your virtual environment
- Try running demos
## Additional Notes
The repository should be updated periodically.
