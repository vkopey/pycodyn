// © Volodymyr Kopei, 2017, email: vkopey@gmail.com
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

model Oscillator // class-model
  Mass mass1(s(start=-1), v(start=0), m=3402.0); // object with initial conditions
  SpringDamper spring1(c=39694.0, d=1856.0); // object
  Fixed fixed1(s0=0); // object
equation // additional equations
  // creates a system of equations (see Flange class)
  connect(fixed1.flange, spring1.flange_a);
  connect(spring1.flange_b, mass1.flange_a);
end Oscillator;
