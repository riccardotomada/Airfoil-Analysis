Airfoil Analysis is a tool which implements the Hess - Smith panel method to compute the aerodynamic coefficients developed by 4 and 5 digits NACA airfoil, in the 2D incompressible steady flow context.

The analysis can be carried out for a single airfoil as well as considering the wing - tail interaction. In both cases, the ground effect condition may be taken into account.

![image1_auto_x2](https://user-images.githubusercontent.com/91890502/189160682-6aae1055-c78a-472e-802d-7b4ab57b54b1.jpg | width=10)
![image2_auto_x2](https://user-images.githubusercontent.com/91890502/189160711-993e2748-de79-4a6e-8c99-66de8b404962.jpg | width=10)
![image3_auto_x2](https://user-images.githubusercontent.com/91890502/189160732-7d5f8acb-8307-458c-ac24-de2dae08966e.jpg | width=10)
![image4_auto_x2](https://user-images.githubusercontent.com/91890502/189160747-eb78492b-b7ed-43d5-b8da-4b94618588b2.jpg | width=10)



Here a list of the main features:

• Geometry definition and visualization

• Aerodynamic coefficients calculation: 
- The pressure coefficient (Cp) is computed pointwise both on the upper and lower surface of the airfoil. Its behaviour is then plotted
- The lift coefficient (Cl) is computed both using the Kutta - Joukowski theorem and the Cp integration procedure
- The pitching moment coefficient (Cm) is computed via Cp integration
- A viscous flow correction method is implemented in order to compute the lift coefficient and the drag coefficient (Cd) in a more accurate way. It consists of the Thwaites' method for laminar flows coupled with the Michel's criterion for transition and the Head's method for the turbulent part of the boundary layer. The skin friction coefficient is retrieved via the Ludwieg-Tillman skin-friction law.

• Thin airfoil theory - rough approximation: get immediate aerodynamic coefficient results without waiting for the panel method output. But remember, it is always an approximation.

• Cl - alpha plot: look how the Cl value changes at different angles of attack. At the moment it is available only for the inviscid case.

• Dataset Generation: create a dataset from a range of parameters defined by you and retrieve their coefficients in a fast and accurate way


This project was born during a fluid dynamics bachelor course laboratory, and it has been revised during the master degree. It is intended for educational purposes only.
