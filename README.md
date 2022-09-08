Airfoil Analysis is a tool which implements the Hess - Smith method to compute the aerodynamic coefficients developed by 4 and 5 digits NACA airfoil, in the 2D incompressible steady flow context.

The analysis can be carried out for a single airfoil as well as considering the wing - tail interaction. In both cases, the ground effect condition may be taken into account.


Here a list of the main features:

• Geometry definition and visualization

• Aerodynamic coefficients calculation: 
- The pressure coefficient (Cp) is computed pointwise both on the upper and lower surface of the airfoil. Its behaviour is then plotted on a graph
- The lift coefficient (Cl) is computed both using the Kutta - Joukowski theorem and the Cp integration way
- The moment coefficient is computed via Cp integration
- A viscous flow correction method is implemented in order to compute the lift coefficient and the drag coefficient (Cd) in a more accurate way (Thwaites method)

• Thin airfoil theory - rough approximation: get immediate aerodynamic coefficient results without waiting for the panel method output. But remember, it is always an approximation.

• Cl - alpha plot: look how the Cl value changes at different angles of attack


This project was born during a fluid dynamics bachelor course laboratory. It is intended for educational purposes.
