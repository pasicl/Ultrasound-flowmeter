Rf = 590;
Rg = 33;
Rs = 55;
R2 = 100;
R1 = 1e3;
In2 = 1.7e-12;
G2 = 2;
bw = 5;

Eni=0.9e-9;
Eni2 = 8e-9;
Ibn = 2e-12;
Ibi=Ibn;
G1 = 1 + Rf/Rg;
Kb = 1.38E-23;
T = 300;

Eo1 = sqrt((Eni*G1)^2+(Ibn*Rs*G1)^2+4*Kb*T*Rs*G1^2+(Ibi*Rf)^2+4*Kb*T*Rf+4*Kb*T*Rf^2/Rg) % Output noise of non-inverting amplifier /sqrt(bw)

Eo2 = sqrt(2*(G2*Eo1)^2 + 4*Kb*T*R2 + Eni2 ^2 * (G2/2)^2 + In2^2*(R2*R1/(R2+R1))^2) %Output noise of line driver

Eadc = 2/2^14/sqrt(12)*sqrt(1/125e6)

Et = sqrt(Eo2^2+Eadc^2)*sqrt(bw)




