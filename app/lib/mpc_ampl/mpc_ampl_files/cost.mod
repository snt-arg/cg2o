# Put these AFTER p01,p02,p10,p11,p20 are set (via let), or define and then "let" them.
param den := (p02*p20 - (p11^2)/4);

param A := (p01*p11)/(4*den);
param B := (p02*p10)/(2*den);
param C := (p01*p20)/(2*den);
param D := (p10*p11)/(4*den);

minimize cost:
  sum {k in K_u} (
     Omega1 * (
          p02 * (f_t[k] + C - D)^2
        + p20 * (v_h[k]^2 - A + B)^2
        + p11 * (v_h[k]^2 - A + B) * (f_t[k] + C - D)
      )
    + BrakePenalty          * f_b[k]^2
    + BrakeChangePenalty    * slack_3[k]^2
    + TractionChangePenalty * slack_2[k]^2
    + Omega6                * slack_1[k]^2
  );



   