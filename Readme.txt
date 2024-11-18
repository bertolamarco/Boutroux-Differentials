The main file is G_functions.m and is the front end for the whole package. 
The syntax is: 
[S,D, e,tR]= G_functions(phi,L,e,Precision, PRINT,  Draw, epsilon, Length )

 Finds the quadratic differential $Q= S^2 Delta/E$ where E = prod (z-e),
 deg S = L, deg Delta = M and \sqrt{Q} \simeq phi' + t_0/z + ..., where
 sqrt{Q} has purely imaginary periods. 
 The quadratic differential is Q = S^2 D/E where E is the poly with roots
 at e and deg S =L, deg D =M, deg E = N, N  = 2L+M+2;
 The argument phi=[t_R, ..., t_1,t_0]; t_0 is the residue of
 \sqrt Q. 
 e=[list of points]
 L= degree of L. If the passed L is too large it is defaulted to L=0, with
 a complaint. 

PRINT=true/false; if true (default false), it will save the picture in the
same directory and also a small .tex wrapper with the same filename

Precision= threshold of the functional  below which we stop the
computation (Default = 1e-8)

Draw, epsilon, Length are parameters for the external function that draws
the critical trajectories. They don't affect the data output. 

OUTPUT:  S= the polynomial S;D = the polynomial Delta; e= the same as the
 one passed; tR= phi(1) (leading coefficient). If needed one can change
 the code to return the polynomials S, D(elta) themselves. 

DEFAULTS: if phi is not defined or empty it defaults to phi=1 (i.e. t_0=1;
Chebotarov case); if e is not defined or empty, it defaults to N=randi(10)+1 random points;
if L is not defined it defaults to L=0; DRAW=true, epsilon = 0.004; Length=20

 Written for Matlab R2024a but possibly working for older versions as
 well. 
 
This package is the accompanying code to https://www.arxiv.org/abs/2408.15234

 August 2024
 (c) Marco Bertola
 email: Marco.Bertola@concordia.ca
