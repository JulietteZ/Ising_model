# Ising model
Monte Carlo simulations for 2D Ising model

* sweep()  
All spins are updated based on Metropolis algorithm at each time stamp.

* afm()  
Change the coordination number to simulate different lattices.  
NOTE: need to update energy() to account for all nearest neighbors when changing coordination number.
