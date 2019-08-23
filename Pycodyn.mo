//Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com
within;
package Pycodyn
extends Modelica.Icons.Package;

connector Flange // class-connector
  Real s; // variable (positions at the flange are equal)
  flow Real f; // variable (sum of forces at the flange is zero)
end Flange;

model Fixed // class-model
  parameter Real s0=0; // parameter (constant in time) 
  Flange flange; // object of class Flange
equation // model equations 
  flange.s = s0; 
end Fixed;

partial model Transl // class-model
  Flange flange_a; // object of class Flange
  Flange flange_b; // object of class Flange
end Transl;

model Mass // class-model
  extends Transl; // inheritance of class Transl
  parameter Real m(min=0, start=1); // parameter 
  Real s; // variable
  Real v(start=0); // variable with initial condition 
  Real a(start=0); // variable with initial condition 
equation // model equations 
  v = der(s);
  a = der(v);
  m*a = flange_a.f + flange_b.f;
  flange_a.s = s;
  flange_b.s = s;
end Mass;

model SpringDamper // class-model
  extends Transl; // inheritance of class Transl
  parameter Real c(final min=0, start=1); // parameter 
  parameter Real d(final min=0, start=1); // parameter 
  Real s_rel(start=0); // variable
  Real v_rel(start=0); // variable
  Real f; // variable
equation // model equations 
  f = c*s_rel+d*v_rel;
  s_rel = flange_b.s - flange_a.s;
  v_rel = der(s_rel);
  flange_b.f = f;
  flange_a.f = -f;
end SpringDamper;

model Motion // class-model
  parameter Real A=2.1/2;
  parameter Real n=6.4/60;
  Flange flange; // object of class Flange
equation // model equations 
  flange.s = A*sin(6.283185307179586*n*time); 
end Motion;

model Force // class-model
  Real f(start=0);
  Flange flange; // object of class Flange
equation // model equations 
  flange.f = f; 
end Force;

model Oscillator // class-model
  Mass mass1(s(start=-1), v(start=0), m=3961.0); // object with initial conditions
  SpringDamper spring1(c=44650.0, d=2120.7); // object
  Fixed fixed1(s0=0); // object
equation // additional equations
  // creates a system of equations (see Flange class)
  connect(fixed1.flange, spring1.flange_a);
  connect(spring1.flange_b, mass1.flange_a);
end Oscillator;

model Pumping // class-model
  Mass mass1(s(start=-0.9), v, m=3961.0); // object with initial conditions
  SpringDamper spring1(c=44650.0, d=2120.7); // object
  Motion motion1(A=2.1/2, n=6.4/60);
  Force force1;
equation // additional equations
  connect(motion1.flange, spring1.flange_a);
  connect(spring1.flange_b, mass1.flange_a);
  connect(mass1.flange_b, force1.flange);
algorithm
  if mass1.v <= 0 then
    force1.f:=34687.0;
  else
    force1.f:=34687.0+18499.0*tanh(abs(mass1.v )/0.01);
  end if;
end Pumping;

end Pycodyn;