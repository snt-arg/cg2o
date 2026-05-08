subject to sub_eq_b  {k in K_u} :
           (1-space_mode) * v_h[k+1] + (space_mode) * v_h[k+1]**2 
           =  
            (1-space_mode) * (
                v_h[k]
              + delta_t/m_eq * (f_t[k] - f_b[k])
              - delta_t/m_eq * (
                    m_total*g*sin(alpha)
                  + c_r*m_total*g*cos(alpha)
                  + c_a * v_h[k]**2
                )
            )
            + space_mode * (
                v_h[k]^2 * (1 - 2*c_a*delta_s/m_eq)
              + 2*delta_s/m_eq * (f_t[k] - f_b[k])
              - 2*delta_s/m_eq * m_total*g * (sin(alpha) + c_r*cos(alpha))
            );

  
                   


subject to sub_eq_c {k in K_u} :
           d_h[k+1]
           = 
           (1-space_mode) * (d_h[k] + 0.5 * delta_t * (v_p[k] + v_p[k+1] - v_h[k] - v_h[k+1]))
            + space_mode * (d_h[k] + ( delta_s / (v_h[k] + v_h[k+1])) * (v_p[k] + v_p[k+1]) - delta_s);
 